// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "DirectCoupledPDE.hh"

#include "BasePairCoupling.hh"
#include "PDE/newmark.hh"

// include solveStep drivers
#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"

// include PDE classes
#include "PDE/SinglePDE.hh"

#include "Domain/domain.hh"
#include "Driver/singleDriver.hh"
#include "Driver/transientdriver.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/resultHandler.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  DirectCoupledPDE::DirectCoupledPDE( Grid *aptgrid, ParamNode* paramNode )

    : StdPDE( aptgrid, paramNode ) {

    ENTER_FCN( "DirectCoupledPDE::DirectCoupledPDE" );

    totalUnknowns_ = 0;
  }


  // **************
  //   Destructor
  // **************
  DirectCoupledPDE::~DirectCoupledPDE() {

    ENTER_FCN( "DirectCoupledPDE::~DirectCoupledPDE" );

    // The following cannot easily be deleted by StdPDE from which we inherit
    // them, since SinglePDE also inherits them.
    delete algsys_;
    delete solveStep_;
    delete TS_alg_;

    // Delete BasePairCoupling objects
    for ( UInt i = 0; i < couplings_.GetSize(); i++ ) {
      delete couplings_[i];
    }
    couplings_.Clear();

    delete solVec_;
    delete rhsVec_;
  }


  // *****************
  //   SetSinglePDEs
  // *****************
  void DirectCoupledPDE::SetSinglePDEs( const StdVector<SinglePDE*> &pdes ) {

    ENTER_FCN( "DirectCoupledPDE::SetSinglePDEs" );

    singlePDEs_ = pdes;
    
    // create pdename
    for ( UInt i = 0; i < singlePDEs_.GetSize()-1; i++ ) {
      pdename_ += singlePDEs_[i]->GetName();
      pdename_ += "-";
    }

    pdename_ += singlePDEs_[singlePDEs_.GetSize()-1]->GetName();
  }


  // ****************
  //   SetCouplings
  // ****************
  void DirectCoupledPDE::SetCouplings( const StdVector<BasePairCoupling*>
                                       &couplings ) {
    ENTER_FCN( "DirectCoupledPDE::SetCouplings" );
    couplings_ = couplings;
  }


  // ********
  //   Init
  // ********
  void DirectCoupledPDE::Init( UInt sequenceStep ) {

    ENTER_FCN( "DirectCoupledPDE::Init" );
  
    sequenceStep_ = sequenceStep;

    // Check, whether we shall generate an SBM_System
    bool genSBMSys = false;
    ParamNode * linSysNode = 
      param->Get( "sequenceStep", "index", GenStr(sequenceStep) )
      ->Get("linearSystems", false );
    if( linSysNode ) {
      ParamNode * specSysNode = linSysNode
        ->Get("system","name","direct", false);
      if( specSysNode ) {
        ParamNode * matrixNode = specSysNode->Get("matrix", false );
        if( matrixNode ) {
          genSBMSys = matrixNode->Has("sbmMatrix");
        }
      }
    }

    // Create algebraic system and pass it to SinglePDEs
    if ( genSBMSys == true ) {
      algsys_ = new SBM_System();
    }
    else {
      algsys_ = new StandardSystem();
    }

    // Get parameter and report object of OLAS
    olasParams_ = algsys_->GetOLASParams();
    olasReport_ = algsys_->GetOLASReport();
  
    // ----------------------------
    //  Detection of analysis type
    // ----------------------------

    analysistype_ = domain->GetSingleDriver()->GetAnalysisType( );

    // Create new Assemble object
    // NOTE: At the moment we can only couple mechanic, piezo and
    // acoustic PDEs directly. All these PDEs have second order
    // time derivatives. This is why we pass a magic '2' directly
    // as maximum time derivative order to assemble.
    // As soon as we have a more sophisticated way to cope with
    // this problem (e.g. register each pde with pdeIdType and
    // maximum time derivative order), we should change this here!
    assemble_ = new Assemble( algsys_, analysistype_, 2 );

    // Initialize timestepping
    if ( analysistype_ == TRANSIENT ) {
      InitTimeStepping();
    }

    // activate direct coupling information
    // and initialize all single pdes
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->algsys_ = algsys_;
      singlePDEs_[i]->assemble_ = assemble_;
      singlePDEs_[i]->SetDirectCoupling();
      // Initialize all SinglePDEs
      singlePDEs_[i]->Init( sequenceStep );
    }

    // Get information about number of dirichlet values,
    // dofs, constraints and needed matrices
  
    // Iterate over all single PDEs and collect data about
    // included boundary conditions
    shared_ptr<EqnMap> eqn;
    for ( UInt i=0; i<singlePDEs_.GetSize(); i++ ) {
      eqn = singlePDEs_[i]->GetEqnMap();
      totalUnknowns_ += eqn->GetNumEqns();
    }
    
    if ( analysistype_ == TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain->GetSingleDriver())
        ->GetDeltaT();
      TS_alg_->Init( dt, totalUnknowns_ );
    }


    // Set correct size of direct solution value
    if ( analysistype_ == HARMONIC  ) {
      solVec_ = new Vector<Complex>;
      rhsVec_ = new Vector<Complex>;
    }
    else {
      solVec_ = new Vector<Double>;
      rhsVec_ = new Vector<Double>;
    }
    
    solVec_->Resize(totalUnknowns_);

    // TEMPORARY CHANGE CHANGE CHANGE
    BaseNodeStoreSol *ptNodeSol;

    if ( analysistype_ == HARMONIC  ) {
      solVec_->Init( Complex(0.0, 0.0 ) );
      Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        singlePDEs_[i]->solVec_ = solVec_;
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }
    } else {
      solVec_->Init( 0.0 );
      Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        singlePDEs_[i]->solVec_ = solVec_;
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }
    }

  
    //! Augment nonlinearity information
    //! from PiezoCoupling, mechPDE and elecPDE

    bool globalNonLin = false;
    bool globalNonLinMaterial = false;
    bool globalNonLinHysteresis = false;

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      if(singlePDEs_[i]->IsNonLin())
        globalNonLin=true;
      if(singlePDEs_[i]->IsNonLinMaterial())
        globalNonLinMaterial=true;
    }
    //    std::cout<< "Direct CoupledPDE - Init globalNonLinMaterial " << globalNonLinMaterial << std::endl;
    
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      //      Is NonLin()
      if(couplings_[i]->nonLin_==true)
        globalNonLin=true;
      if(couplings_[i]->nonLinMaterial_==true)
        globalNonLinMaterial=true;

      if(couplings_[i]->nonLinHysteresis_==true) {
        globalNonLinHysteresis = true;
        isHysteresis_ = true;
      }
    }
    
    nonLin_=globalNonLin;    
    nonLinMaterial_=globalNonLinMaterial;    

    if ( !globalNonLinHysteresis ) {
      // copy nonlinearity information to singlePDEs
      for (UInt i=0; i<singlePDEs_.GetSize(); i++){
        singlePDEs_[i]->SetNonLinearity(globalNonLin);
        singlePDEs_[i]->SetMaterialNonLinearity(globalNonLin);
      }
      
      // copy nonlinearity information to couplings   
      for (UInt i=0; i<couplings_.GetSize(); i++) {
        couplings_[i]->SetNonLinearity(globalNonLin);    
        couplings_[i]->SetMaterialNonLinearity(globalNonLin);
      }
    }
    
    // define solveStep-driver
    DefineSolveStep();

    // Pass SolveStep object to all single pdes
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->solveStep_ = solveStep_;
    }
    // Initialize all Coupling Objects 
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->SetAlgSys( algsys_ );
      couplings_[i]->SetAssemble( assemble_ );
      couplings_[i]->Init( sequenceStep_ );
    }

  }


  // **************************
  //   WriteGeneralPDEdefines
  // **************************
  void DirectCoupledPDE::WriteGeneralPDEdefines() {
    ENTER_FCN( "DirectCoupledPDE::WriteGeneralPDEdefines" );

    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->WriteGeneralPDEdefines();
    }
  }


  // **********
  //   SetBCs
  // **********
  void DirectCoupledPDE::SetBCs( const Double atimestep ) {

    ENTER_FCN( "DirectCoupledPDE::SetBCs" );

    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->SetBCs( atimestep );
    }
  }



  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double DirectCoupledPDE::RhsL2Norm(Vector<Double>& actRHS)
  {
    ENTER_FCN( "DirectCoupledPDE::RhsL2Norm" );

    Integer eqnNr;

    for ( UInt j = 0; j < singlePDEs_.GetSize(); j++ ) {

      IdBcList idbcList = singlePDEs_[j]->GetIDBCList();

      // Eliminate inhom. dirichlet node from RHS (due to penalty formulation)
      for ( UInt i = 0; i < idbcList.GetSize(); i++ ) {
        
        // Get grip of current idBC
        InhomDirichletBc const & actBc = *idbcList[i];

        // Get entity iterator
        EntityIterator it = actBc.entities->GetIterator();

        
        for ( it.Begin(); !it.IsEnd(); it++ ) {
          eqnNr = singlePDEs_[j]->GetEqnMap()->GetEqn( *actBc.result, it, actBc.dof );
          if ( eqnNr != 0 ) {
            actRHS[(eqnNr-1)] = 0.0;
          }
        }
      } 
    }
    
    return actRHS.NormL2();
  }


  // ****************
  //   DefineAlgSys
  // ****************
  void DirectCoupledPDE::DefineAlgSys() {

    ENTER_FCN( "DirectCoupledPDE::DefineAlgSys" );
 
    
    std::string pdeName;
    PdeIdType pdeId;
    shared_ptr<EqnMap> eqn;
    
    // Set linear system parameters for OLAS
    //
    // NOTE: Using current naming conventions in the XML Schema definitions
    //       the linear system for a direct coupled PDE problem is always
    //       called "direct", thus there can currently only be one of them.
    ReadOlasParams( "direct" );

    // Begin setup of the matrix graph
    algsys_->GraphSetupInit( singlePDEs_.GetSize() );

    // iterate over all singlePDE and register them
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {

      // obtain PDE identification tag from algebraic system
      // and set number of dirichlet and constraint equations
      pdeName= singlePDEs_[i]->GetName();
      eqn = singlePDEs_[i]->GetEqnMap();
      pdeId = singlePDEs_[i]->GetPDEId();
      algsys_->RegisterPDE( pdeId, eqn->GetNumEqns(),
                            eqn->GetNumLastFreeDof() );

      // Let the PDE set its Dirichlet information and related stuff
      singlePDEs_[i]->DefineAlgSys();
    }

    // CURRENTLY NOT REQUIRED
    //
    // // iterate over all coupling objects and register them
    // for ( UInt i = 0; i < couplings_.GetSize(); i++ ) {
    // 
    //   // register forward coupling
    //   algsys_->RegisterCoupling( couplings_[i]->GetPdeId1(),
    //                              couplings_[i]->GetPdeId2() );
    // 
    //   // register backward coupling
    //   algsys_->RegisterCoupling( couplings_[i]->GetPdeId2(),
    //                              couplings_[i]->GetPdeId1() );
    // }

    // iterate over all singlePDE and setup matrix graph
    // trigger the creation and assembly of the matrix graph
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      PdeIdType id = singlePDEs_[i]->GetPDEId();
      algsys_->AssembleInit( id, id, false );
      assemble_->SetupMatrixGraph( id, id );
      algsys_->AssembleDone( id, id, false );
    }

    // For SBM_Systems we must obtain the re-orderings now
    // and pass it to the EQN-objects
    if ( dynamic_cast<SBM_System*>(algsys_) != NULL ) {
      IncorporateReordering();
    }

    // Setup matrix graph of coupling objects
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      PdeIdType id1 = couplings_[i]->GetPdeId1();     
      PdeIdType id2 = couplings_[i]->GetPdeId2();

      // setup matrix graph
      algsys_->AssembleInit( id1, id2, true );
      assemble_->SetupMatrixGraph( id1, id2 );
      algsys_->AssembleDone( id1, id2, true );

      // Note: also the second part has to be assembled!
      algsys_->AssembleInit( id2, id1, true );
      assemble_->SetupMatrixGraph( id2, id1 );
      algsys_->AssembleDone( id2, id1, true );
    }

    // Finish assembly of the matrix graph
    algsys_->GraphSetupDone();

    // For StandardSystems we must obtain the re-orderings at this point
    // and pass it to the EQN-objects
    if ( dynamic_cast<StandardSystem*>(algsys_) != NULL ) {
      IncorporateReordering();
    }

#ifdef DEBUG
    // Print information from assemble class
    assemble_->PrintInfo( *debug );
#endif

    // Allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();

    for( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      if (singlePDEs_[i]->memento_ != NULL &&
          singlePDEs_[i]->mementoUsage_ == PDEMemento::START_VALUE ) {
        singlePDEs_[i]->IncorporateMemento();
      }
    }
        
  }

  void DirectCoupledPDE::SaveSolution( const Double * ptSol, UInt size) {
    ENTER_FCN( "DirectCoupledPDE::SaveSolution" );

    BaseNodeStoreSol *ptNodeSol;
    Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVec_);
    solHelp.Resize(size);
    
    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }
    
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      // set pointer to solution object of the PDE
      ptNodeSol = singlePDEs_[i]->getPDESolution();
      ptNodeSol->SetAlgSysDataPointer( size, solHelp.GetPointer() );
      singlePDEs_[i]->solVec_  = solVec_;
      
    }
  }

  void DirectCoupledPDE::SaveSolution( const Complex * ptSol, UInt size) {
    ENTER_FCN( "DirectCoupledPDE::SaveSolution" );

    BaseNodeStoreSol *ptNodeSol;
    Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*solVec_);
    solHelp.Resize(size);
    
    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }
    
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      // set pointer to solution object of the PDE
      ptNodeSol = singlePDEs_[i]->getPDESolution();
      ptNodeSol->SetAlgSysDataPointer( size, solHelp.GetPointer() );
      
    }
  }

   void DirectCoupledPDE::SaveRHS( const Double * ptSol, UInt size) {
    ENTER_FCN( "DirectCoupledPDE::SaveRhs" );

    Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*rhsVec_);
    solHelp.Resize(size);
    
    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }
    
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      // set pointer to solution object of the PDE
      singlePDEs_[i]->rhsVec_  = rhsVec_;
      
    }
  }
  void DirectCoupledPDE::SaveRHS( const Complex * ptSol, UInt size) {
    ENTER_FCN( "DirectCoupledPDE::SaveRhs" );

    Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*rhsVec_);
    solHelp.Resize(size);
    
    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }
    
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      // set pointer to solution object of the PDE
      singlePDEs_[i]->rhsVec_  = rhsVec_;
      
    }
  }
  

  // ********************
  //   InitTimeStepping
  // ********************
  void DirectCoupledPDE::InitTimeStepping() {

    ENTER_FCN( "DirecCoupledPDE::InitTimeStepping" );

    // Hard Coded
    TS_alg_ = new Newmark( algsys_ );

    // Pass time stepping object to single pdes
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->TS_alg_ = TS_alg_;
    }
  }


  // ======================================================
  // POSTPROC SECTION
  // ======================================================


  void DirectCoupledPDE::WriteRestart( ) 
  {
    ENTER_FCN( "DirectCoupledPDE::WriteRestart" );

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->WriteRestart( );
    }
  }

  void DirectCoupledPDE::ReadRestart( UInt &startStep ) 
  {
    ENTER_FCN( "DirectCoupledPDE::ReadRestart" );

    StdVector<UInt> startSteps( singlePDEs_.GetSize() );
     
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->ReadRestart(startSteps[i]);
    }
     
    for( UInt i = 1; i < startSteps.GetSize(); i++ ) {
      if( startSteps[i] != startSteps[0] ) {
        std::stringstream errMsg;
        errMsg <<  "Error during read in of restart files:\n"
               << "Restart step numbers differ for the different PDEs!\n";
        for( UInt j = 0; j< singlePDEs_.GetSize(); j++ ) {
          errMsg << singlePDEs_[i]->GetName() << "\t"
                 << startSteps[i] << "\n";
        }
        EXCEPTION( errMsg.str().c_str() );
      }
       
    }
  }
  

  void DirectCoupledPDE::WriteResultsInFile(const UInt kstep,
                                            const Double asteptime,
                                            UInt stepOffset,
                                            Double timeOffset)
  {
    ENTER_FCN( "DirectCoupledPDE::WriteResultsInFile" );

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->WriteResultsInFile( kstep, asteptime, 
                                          stepOffset, timeOffset);
    }
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->WriteResultsInFile( kstep, asteptime, 
                                         stepOffset, timeOffset);
    }

  }

  void DirectCoupledPDE::Finalize() {
    ENTER_FCN( "DirectCoupledPDE::Finalize" );

    for ( UInt i=0; i<singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->Finalize();
    }
    
    for( UInt i=0; i<couplings_.GetSize(); i++ ) {
      couplings_[i]->Finalize();
    }
    

    domain->GetResultHandler()->Finalize();
  }
  
  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void DirectCoupledPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "DirectCoupledPDE::InitCoupling" );
    isIterCoupled_ = true;
  }
  void DirectCoupledPDE::ResetCoupling()
  {
    ENTER_FCN( "DirectCoupledPDE::ResetCoupling ");

    iterCoupledCounter_ = 0;
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->ResetCoupling();
      }
  }
           
  void DirectCoupledPDE::CalcInputCoupling()
  {
    ENTER_FCN( "DirectCoupledPDE::CalcInputCoupling" );
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->CalcInputCoupling();
      }
  }

  void DirectCoupledPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "DirectCoupledPDE::CalcOutputCoupling" );
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->CalcOutputCoupling();
      }
  }




  void DirectCoupledPDE::DefineSolveStep() {
    ENTER_FCN( "DirectCoupledPDE::DefineSolveStep" );
  
    solveStep_ = new StdSolveStep(*this);
  }

  // *************************
  //   IncorporateReordering
  // *************************
  void DirectCoupledPDE::IncorporateReordering() {

    ENTER_FCN( "DirectCoupledPDE::IncorporateReordering" );

    PdeIdType pdeId   = NO_PDE_ID;
    shared_ptr<EqnMap> eqn;
    Integer *newOrder = NULL;

    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      pdeId    = singlePDEs_[i]->GetPDEId();
      eqn      = singlePDEs_[i]->GetEqnMap();
      newOrder = algsys_->GetReordering( pdeId );
      if( newOrder == NULL ) {
        std::cerr << "performing no reordering!";
      }
      
      eqn->ReorderMapping( &newOrder );
    }
  }

}
