#include "DirectCoupledPDE.hh"

#include "BasePairCoupling.hh"
#include "PDE/newmark.hh"

// include solveStep drivers
#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"

// include PDE classes
#include "PDE/SinglePDE.hh"

#include "Domain/domain.hh"
#include "Driver/basedriver.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  DirectCoupledPDE::DirectCoupledPDE( Grid *aptgrid, WriteResults *aOutFile, 
                                      TimeFunc *aTimeFunc )

    : StdPDE( aptgrid, aOutFile, aTimeFunc ) {

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
  void DirectCoupledPDE::Init( UInt bcSequenceStep,
                               std::string bcSequenceTag ) {

    ENTER_FCN( "DirectCoupledPDE::Init" );
  
    bcSequenceIndex_ = bcSequenceStep;
    bcSequenceTag_ = bcSequenceTag;

    // Check, whether we shall generate an SBM_System
    StdVector<std::string> aux;
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    bool genSBMSys = false;
    keyVec  = "linearSystems", "system", "sbmMatrix", "entry";
    attrVec = "", "name", "";
    valVec  = "", "direct", "";
    params->GetList( keyVec, attrVec, valVec, aux );
    genSBMSys = aux.GetSize() == 1;

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

    analysistype_ = domain->GetDriver()->GetAnalysisType( pdename_ );

    // Create new Assemble object
    assemble_ = new Assemble( algsys_, analysistype_, ptTimeFunc_ );

    // activate direct coupling information
    // and initialize all single pdes
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->algsys_ = algsys_;
      singlePDEs_[i]->assemble_ = assemble_;
      singlePDEs_[i]->SetDirectCoupling();
      // Initialize all SinglePDEs
      singlePDEs_[i]->Init( bcSequenceStep, bcSequenceTag );

      
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
    
    // Initialize timestepping
    if ( analysistype_ == TRANSIENT ) {
      InitTimeStepping();
    }


    // Set correct size of direct solution value
    if ( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ) {
      solVec_ = new Vector<Complex>;
    }
    else {
      solVec_ = new Vector<Double>;
    }
    
    solVec_->Resize(totalUnknowns_);


    // TEMPORARY CHANGE CHANGE CHANGE
    BaseNodeStoreSol *ptNodeSol;

    if ( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ) {
      solVec_->Init( Complex(0.0, 0.0 ) );
      Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }
    } else {
      solVec_->Init( 0.0 );
      Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }
    }

    // define solveStep-driver
    DefineSolveStep();
    
    // Pass SolveStep object to all single pdes
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->solveStep_ = solveStep_;
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

    // Initialize all Coupling Objects and setup their matrix graph
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->SetAlgSys( algsys_ );
      couplings_[i]->SetAssemble( assemble_ );
      couplings_[i]->Init( bcSequenceIndex_, bcSequenceTag_ );
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

  void DirectCoupledPDE::PostProcess()
  {
    ENTER_FCN( "DirectCoupledPDE::PostProcess" );
  
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->PostProcess();
    }

    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->PostProcess();
    }
  }

   void DirectCoupledPDE::WriteRestart(const UInt nstep, UInt totalUnknowns) 
   {
     ENTER_FCN( "DirectCoupledPDE::WriteRestart" );

     singlePDEs_[0]->WriteRestart(nstep,totalUnknowns_);
   }

   void DirectCoupledPDE::ReadRestart(UInt &startStep, UInt totalUnknowns) 
   {
     ENTER_FCN( "DirectCoupledPDE::ReadRestart" );

     singlePDEs_[0]->ReadRestart(startStep,totalUnknowns_);
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

  void DirectCoupledPDE::WriteHistoryInFile(const UInt kstep,
                                            const Double asteptime,
                                            UInt stepOffset,
                                            Double timeOffset)
  {
    ENTER_FCN( "DirectCoupledPDE::WriteHistoryInFile" );

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->WriteHistoryInFile( kstep, asteptime, 
                                          stepOffset, timeOffset);
    }

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
