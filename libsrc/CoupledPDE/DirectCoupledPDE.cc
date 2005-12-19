#include "DirectCoupledPDE.hh"

#include "BasePairCoupling.hh"
#include "PDE/newmark.hh"

// include solveStep drivers
#include "Driver/stdSolveStep.hh"
#include "Driver/solveStepMech.hh"
#include "Driver/solveStepAcousticMechBubble.hh"
// include PDE classes
#include "PDE/SinglePDE.hh"

#include "PDE/nodeEQN.hh"
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
  
    // activate direct coupling information
    // and initialize all single pdes
    std::set<AnalysisType> analysisTypes;
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->SetAlgebraicSystem(algsys_);
      singlePDEs_[i]->SetDirectCoupling();
    
      // Initialize all SinglePDEs
      singlePDEs_[i]->Init( bcSequenceStep, bcSequenceTag );

      // Insert analysistype into set
      analysisTypes.insert(singlePDEs_[i]->GetAnalysisType());
    }

    // ----------------------------
    //  Detection of analysis type
    // ----------------------------

    // NOTE: The current implementation is not sophisticated. We need
    //       to develop a better concept here
    if ( analysisTypes.find(TRANSIENT) != analysisTypes.end() ) {
      analysistype_ = TRANSIENT;
    }
    else if ( analysisTypes.find(STATIC) != analysisTypes.end() ) {
      analysistype_ = STATIC;
    }
    else if ( analysisTypes.find(HARMONIC) != analysisTypes.end() ) {
      analysistype_ = HARMONIC;
    }
    else {
      (*error) << "I was not able to determine the correct analysistype for "
               << "your set of PDEs!";
      Error( __FILE__, __LINE__ );
    }


    // We simply take the assemble object of the first SinglePDE
    assemble_ = singlePDEs_[0]->getPDE_assemble();

    // Get information about number of dirichlet values,
    // dofs, constraints and needed matrices
  
    // Iterate over all single PDEs and collect data about
    // included boundary conditions
    UInt numDroppedDofs = 0;
    NodeEQN * eqn;
    for ( UInt i=0; i<singlePDEs_.GetSize(); i++ ) {
      eqn = singlePDEs_[i]->getPDE_eqnData();
      numDroppedDofs += eqn->GetNumDroppedDofs();
      totalUnknowns_ += eqn->GetNumEQNs() * eqn->GetNumDofsPerEQN();
    }
    
    // Initialize timestepping
    if ( analysistype_ == TRANSIENT ) {
      InitTimeStepping();
    }

    // define solveStep-driver
    DefineSolveStep();
    
    // Pass SolveStep object to all single pdes
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->SetSolveStep( solveStep_ );
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
      Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }
    } else {
      Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }
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
    NodeEQN *eqn = NULL;

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
      eqn = singlePDEs_[i]->getPDE_eqnData();
      pdeId = singlePDEs_[i]->GetPDEId();
      algsys_->RegisterPDE( pdeId, eqn->GetNumEQNs(),
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
      singlePDEs_[i]->SetupMatrixGraph();
    }

    // For SBM_Systems we must obtain the re-orderings now
    // and pass it to the EQN-objects
    if ( dynamic_cast<SBM_System*>(algsys_) != NULL ) {
      IncorporateReordering();
    }

    // Initialize all Coupling Objects and setup their matrix graph
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->SetAlgSys( algsys_ );
      couplings_[i]->Init( bcSequenceIndex_, bcSequenceTag_ );
      couplings_[i]->SetupMatrixGraph();
    }

    // Finish assembly of the matrix graph
    algsys_->GraphSetupDone();
  
    // For StandardSystems we must obtain the re-orderings at this point
    // and pass it to the EQN-objects
    if ( dynamic_cast<StandardSystem*>(algsys_) != NULL ) {
      IncorporateReordering();
    }

    // Allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();
  }


  UInt DirectCoupledPDE::GetNumRestraints() {
    ENTER_FCN( "DirectCoupledPDE::GetNumRestraints" );

    UInt totalNumRestraints = 0;
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      totalNumRestraints += singlePDEs_[i]->GetNumRestraints();
  
    return totalNumRestraints;
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
    TS_alg_ = new Newmark( algsys_, totalUnknowns_ );
    // Pass time stepping object to single pdes
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->SetTimeStepping(TS_alg_);
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
    isIterCoupled_ = TRUE;
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


  // ======================================================
  // METHODS FOR ASSEMBLING
  // ======================================================


  void DirectCoupledPDE::AssembleMatrices() {

    // Assembly of diagonal-matrices
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->AssembleMatrices();
    }

    // Assembly of off-diagonal entries (coupling objcts)
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->AssembleMatrices();
    }
  }

  void DirectCoupledPDE::AssembleSrcRHS(const Double time) {
    ENTER_FCN( "DirectCoupledPDE::AssembleSrcRHS" );

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->AssembleSrcRHS( time );
      }
  
  }

  void DirectCoupledPDE::AssembleNLRHS(const Double time) {
    ENTER_FCN( "DirectCoupledPDE::AssembleNLRHS" );

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->AssembleNLRHS( time );
      }
  }

  void DirectCoupledPDE::AssembleSprings( const Double time) {
    ENTER_FCN( "DirectCoupledPDE::AssembleSprings" );

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->AssembleSprings(time );
      }
  }

  void DirectCoupledPDE::AssembleSpecial( ) {
    ENTER_FCN( "DirectCoupledPDE::AssembleSpecial" );
    
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->AssembleSpecial( );
      }
  }
  
  void DirectCoupledPDE::InitNonLinMatrices() {
    ENTER_FCN( "DirectCoupledPDE::InitNonLinMatrices" );

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) 
      {
        singlePDEs_[i]->InitNonLinMatrices();
      }
  }

  void DirectCoupledPDE::SetReassemble() {

    ENTER_FCN( "DirectCoupledPDE::SetReassemble" );

    for (UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->SetReassemble();
    }
    for (UInt i = 0; i < couplings_.GetSize(); i++ ) {
      couplings_[i]->SetReassemble();
    }
  }

 void DirectCoupledPDE::SetFrequency(Double actFreq) {

    ENTER_FCN( "DirectCoupledPDE::SetFrequency" );
    
    for (UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->SetFrequency(actFreq);
    }
    
    for (UInt i = 0; i < couplings_.GetSize(); i++ ) {
      couplings_[i]->SetFrequency(actFreq);
    }
 }
  
  void DirectCoupledPDE::SetupMatrixGraph() {
    // nothing to do here, since this method gets only called for 
    // SinglePDEs
  }

  void DirectCoupledPDE::DefineSolveStep() {
    ENTER_FCN( "DirectCoupledPDE::DefineSolveStep" );
  
    // HARD CODED
    // In case of acou-mech-coupling with pressure formulation
    // use as driver SolveStepAcousticMechBubble
    // otherwise SolveStepMech
    bool acouMechBubble = false;
    std::string pdeName;
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> auxVec;
    
    for (UInt i = 0 ; i < singlePDEs_.GetSize() ; i++){
      pdeName = singlePDEs_[i]->GetName();  
      if ( pdeName == "acoustic" ){
	
	BubbleDynType bubbleDynType = NOBUBBLETYPE;
	keyVec = "acoustic", "bubbles", "bubbleType";
	attrVec = "", "tag";
	valVec =  "", bcSequenceTag_;
	params->GetList(keyVec, attrVec, valVec, auxVec);
	if ( auxVec.GetSize() == 1 ) {
	  String2Enum( auxVec[0], bubbleDynType );
	}
	else {
	  bubbleDynType = NOBUBBLETYPE;
	}
	if ( bubbleDynType !=  NOBUBBLETYPE){
	  std::string str;
	  SolutionType formulation;
	  params->Get( "formulation", str, "pdeList", "acoustic" );
	  String2Enum( str, formulation );
	  if ( formulation ==  ACOU_PRESSURE){
	    acouMechBubble = true;
	  }
	}
	break;
      }
    }
    if ( acouMechBubble == true ){
      solveStep_ = new SolveStepAcousticMechBubble(*this, GILMORE );
    }
    else{
      solveStep_ = new SolveStepMech(*this);
    }
  }

  // *************************
  //   IncorporateReordering
  // *************************
  void DirectCoupledPDE::IncorporateReordering() {

    ENTER_FCN( "DirectCoupledPDE::IncorporateReordering" );

    PdeIdType pdeId   = NO_PDE_ID;
    NodeEQN *eqn      = NULL;
    Integer *newOrder = NULL;

    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      pdeId    = singlePDEs_[i]->GetPDEId();
      eqn      = singlePDEs_[i]->getPDE_eqnData();
      newOrder = algsys_->GetReordering( pdeId );
      eqn->ReorderMapping( &newOrder );
    }
  }

}
