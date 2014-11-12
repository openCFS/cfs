// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <sstream>
#include <string>

#include "BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DirectCoupledPDE.hh"
#include "Domain/bcs.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "Driver/assemble.hh"
#include "Driver/singleDriver.hh"
#include "Driver/solveStepAcoustic.hh"
#include "Driver/solveStepMicroPiezo.hh"
#include "Driver/solveStepPiezo.hh"
// include solveStep drivers
#include "Driver/stdSolveStep.hh"
#include "Driver/transientdriver.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/vector.hh"
#include "OLAS/algsys/basesystem.hh"
#include "OLAS/algsys/sbmsystem.hh"
#include "OLAS/algsys/standardsys.hh"
// include PDE classes
#include "PDE/SinglePDE.hh"
#include "PDE/basePDE.hh"
#include "PDE/eqnMap.hh"
#include "PDE/newmark.hh"
#include "PDE/timestepping.hh"
#include "Utils/basenodestoresol.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
class Grid;
class PDECoupling;

  DirectCoupledPDE::DirectCoupledPDE( Grid *aptgrid, PtrParamNode paramNode )

    : StdPDE( aptgrid, paramNode ) {


    totalUnknowns_ = 0;
  }


  // **************
  //   Destructor
  // **************
  DirectCoupledPDE::~DirectCoupledPDE() {

    // The following cannot easily be deleted by StdPDE from which we inherit
    // them, since SinglePDE also inherits them.
    if( algsys_ ) {
      delete algsys_;
      algsys_ = NULL;
    }
    if( solveStep_ ) {
      delete solveStep_;
      solveStep_ = NULL;
    }
    if( TS_alg_ ){
      delete TS_alg_;
      TS_alg_ = NULL;
    }

    // Delete BasePairCoupling objects
    for ( UInt i = 0; i < couplings_.GetSize(); i++ ) {
      delete couplings_[i];
    }
    couplings_.Clear();

    if( solVec_ ) {
      delete solVec_;
      solVec_ = NULL;
    }
    if( rhsVec_ ) {
      delete rhsVec_;
      rhsVec_ = NULL;
    }
    if ( needSolPrev_ ) {
      delete solVecPrev_;
      solVecPrev_ = NULL;
    }
    
    
    if(assemble_) {
      delete assemble_;
      assemble_ = NULL;
    }
  }


  // *****************
  //   SetSinglePDEs
  // *****************
  void DirectCoupledPDE::SetSinglePDEs( const StdVector<SinglePDE*> &pdes ) {

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
    couplings_ = couplings;
  }

  // ****************
  //   SetInitial conditions
  //   from the pde members
  // ****************

  void DirectCoupledPDE::SetInitialCondition() {

    Vector< Double > aux;
    aux.Init();

    // Construct the initial solution vector
    shared_ptr<EqnMap> eqn;
    UInt singleUnknowns=0;
    UInt lastIndex=0;


    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {

      //------------------------------------------------
      // get the id of the this pde
      //pdeId = singlePDEs_[i]->GetPDEId();
      //------------------------------------------------

      // the PDEs members haven't set up their initial conditions yet
      // -> set the initial conditions of this pde
      singlePDEs_[i]->SetInitialCondition();
      if(singlePDEs_[i]->IsSetInitialCondition()==true){
        this->isSetInitialCondition_=true;
      }


      //------------------------------------------------
      // get the number of unknowns of this pde
      eqn = singlePDEs_[i]->GetEqnMap();
      singleUnknowns = eqn->GetNumEqns();

      // check setup of linear system
      if(singlePDEs_[i]->usePenalty_==false){
        // uses elimination of Inhomogeneous DBC
        //std::cout << "Num of Inhomogeneous DBC = "<< eqn->GetNumInHomDirichletEqns () << std::endl;
        singleUnknowns -= eqn->GetNumInHomDirichletEqns();
        singleUnknowns -= eqn->GetNumInHomDirichletFileEqns();
      }
      //------------------------------------------------

      // init the aux vector
      // **** it is supoussed that i=pdeID ****
      // -> the order of the solution vector is the same as ordering of pdeID
      for (UInt ii = lastIndex; ii < lastIndex+singleUnknowns; ii++) {
        //aux[ii]=singlePDEs_[i]->getInitialCondition();
        aux.Push_back(singlePDEs_[i]->getInitialCondition());
      }

      lastIndex+=singleUnknowns;

    }

    // now we have our solution vector initialized
    //std::cout << "\n al final aux = "<< aux.Serialize() << std::endl;

    if(this->IsSetInitialCondition()==true){


      // save the initial solution vector into algsys
      algsys_->InitSol(aux);

      // save the initial vector in each pde solution vector
      SaveSolution(aux.GetPointer(), aux.GetSize());

      // save the initial solution vector into algsys
      algsys_->InitSol( aux );

    }

  }




  // ********
  //   Init
  // ********
  void DirectCoupledPDE::Init( UInt sequenceStep )
  {
    sequenceStep_ = sequenceStep;

    infoNode_ = info->Get("PDE")->Get("directCoupledPDE", ParamNode::APPEND);
    infoNode_->Get(ParamNode::HEADER)->Get("sequeceStep")->SetValue(sequenceStep);

    // Check, whether we shall generate an SBM_System
    bool genSBMSys = false;
    PtrParamNode linSysNode =
      param->GetByVal( "sequenceStep", std::string("index"), sequenceStep)
      ->Get("linearSystems", ParamNode::PASS );
    if( linSysNode ) {
      PtrParamNode specSysNode = linSysNode
        ->GetByVal("system","name","direct", ParamNode::PASS);
      if( specSysNode ) {
        if( specSysNode->Has("sbmMatrix") ) {
          genSBMSys = true;
        }

      }
    }

    // Create algebraic system and pass it to SinglePDEs
    if ( genSBMSys == true ) {
      algsys_ = new SBM_System(FindLinearSystem("direct"));
    }
    else {
      algsys_ = new StandardSystem(FindLinearSystem("direct"));
    }

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
    if ( analysistype_ == BasePDE::TRANSIENT ) {
      InitTimeStepping();
    }

    // activate direct coupling information
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->SetDirectCoupling();
    }
    
    // activate direct coupling information
    // and initialize all single pdes
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->algsys_ = algsys_;
      singlePDEs_[i]->assemble_ = assemble_;
      // Initialize all SinglePDEs
      singlePDEs_[i]->Init( sequenceStep, infoNode_);

      // check if single PDE really needs previous solution
      if ( singlePDEs_[i]->BelongsPDE2PiezoHyst() 
           || singlePDEs_[i]->BelongsPDE2MicroPiezo() ) {
        needSolPrev_ = true;
      }
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

    if ( analysistype_ == BasePDE::TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain->GetSingleDriver())
        ->GetDeltaT();
      TS_alg_->Init( dt, totalUnknowns_ );
    }


    // Set correct size of direct solution value
    if ( analysistype_ == BasePDE::HARMONIC  ) {
      solVec_ = new Vector<Complex>;
      rhsVec_ = new Vector<Complex>;
    }
    else {
      solVec_ = new Vector<Double>;
      rhsVec_ = new Vector<Double>;
      if ( needSolPrev_ ) {
        solVecPrev_ = new Vector<Double>;
        solVecPrev_->Resize(totalUnknowns_);
      }
    }

    solVec_->Resize(totalUnknowns_);

    // TEMPORARY CHANGE CHANGE CHANGE
    BaseNodeStoreSol *ptNodeSol, *ptNodeSolPrev;

    if ( analysistype_ == BasePDE::HARMONIC  ) {
      solVec_->Init( );
      Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        singlePDEs_[i]->solVec_ = solVec_;
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }
    } else {
      solVec_->Init( );
      Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVec_);
      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        singlePDEs_[i]->solVec_ = solVec_;
        ptNodeSol = singlePDEs_[i]->getPDESolution();
        ptNodeSol->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
      }

      if (  needSolPrev_ ) {
        solVecPrev_->Init( );
        Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVecPrev_);
        for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
          singlePDEs_[i]->solVecPrev_ = solVecPrev_;
          ptNodeSolPrev = singlePDEs_[i]->getPDESolutionPrev();
          ptNodeSolPrev->SetAlgSysDataPointer(totalUnknowns_, solHelp.GetPointer());
        }
      }
    }



    //! Augment nonlinearity information
    //! from PiezoCoupling, mechPDE and elecPDE

    bool globalNonLin = false;
    bool globalNonLinHysteresis = false;

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      if(singlePDEs_[i]->IsNonLin())
        globalNonLin=true;
    }

    for (UInt i=0; i<couplings_.GetSize(); i++) {
      //      Is NonLin()
      if(couplings_[i]->nonLin_==true)
        globalNonLin=true;

      if(couplings_[i]->nonLinHysteresis_==true) {
        globalNonLinHysteresis = true;
        isHysteresis_ = true;
      }
    }

    nonLin_=globalNonLin;

    if ( !globalNonLinHysteresis ) {
      // copy nonlinearity information to singlePDEs
      for (UInt i=0; i<singlePDEs_.GetSize(); i++){
        singlePDEs_[i]->SetNonLinearity(globalNonLin);
      }

      // copy nonlinearity information to couplings
      for (UInt i=0; i<couplings_.GetSize(); i++) {
        couplings_[i]->SetNonLinearity(globalNonLin);
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

    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->WriteGeneralPDEdefines();
    }
  }


  // **********
  //   SetBCs
  // **********
  void DirectCoupledPDE::SetBCs() {


    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      singlePDEs_[i]->SetBCs();
    }
  }



  // calculates L2-norm of RHS regarding dirichlet entries due to penalty
  // formulation by setting them 0
  Double DirectCoupledPDE::RhsL2Norm(Vector<Double>& actRHS)
  {

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



    std::string pdeName;
    PdeIdType pdeId;
    shared_ptr<EqnMap> eqn;

    // Set linear system parameters for OLAS
    //
    // NOTE: Using current naming conventions in the XML Schema definitions
    //       the linear system for a direct coupled PDE problem is always
    //       called "direct", thus there can currently only be one of them.
    ReadOlasParams( "direct" );
    olasInfo_ = info->Get("OLAS")->Get("direct");

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

      // setup matrix graph for upper diagonal(s)
      algsys_->AssembleInit( id1, id2, false );
      assemble_->SetupMatrixGraph( id1, id2 );
      algsys_->AssembleDone( id1, id2, false );

      // setup matrix graph for lower diagonal(s)
      algsys_->AssembleInit( id2, id1, false );
      assemble_->SetupMatrixGraph( id2, id1 );
      algsys_->AssembleDone( id2, id1, false );
    }

    // Finish assembly of the matrix graph
    algsys_->GraphSetupDone();

    // For StandardSystems we must obtain the re-orderings at this point
    // and pass it to the EQN-objects
    if ( dynamic_cast<StandardSystem*>(algsys_) != NULL ) {
      IncorporateReordering();
    }

    // Print information from assemble class
    assemble_->ToInfo(infoNode_->Get(ParamNode::HEADER)->Get("integrators"));

    // Allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();

    // =====================================================================
    // Set the initial conditions
    // =====================================================================
    if ( analysistype_ == TRANSIENT ){
      SetInitialCondition();
    }


    for( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      if (singlePDEs_[i]->memento_ != NULL &&
          singlePDEs_[i]->mementoAsDirichlet_ == false ) {
        singlePDEs_[i]->IncorporateMemento();
      }
    }

  }

  SingleVector* DirectCoupledPDE::GetSolutionVector() {
    return singlePDEs_[0]->GetSolutionVector();
  }

  SingleVector* DirectCoupledPDE::GetPrevSolutionVector() {
    return singlePDEs_[0]->GetPrevSolutionVector();
  }

  void DirectCoupledPDE::SaveSolution( const Double * ptSol, UInt size) {

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
  void DirectCoupledPDE::SavePrevSolution( const Double * ptSolPrev, UInt size) {

    if ( needSolPrev_ ) {
      BaseNodeStoreSol *ptNodeSol;
      Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVecPrev_);
      solHelp.Resize(size);

      for ( UInt i = 0; i < size; i++ ) {
        solHelp[i] = ptSolPrev[i];
      }

      for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
        // set pointer to solution object of the PDE
        ptNodeSol = singlePDEs_[i]->getPDESolutionPrev();
        ptNodeSol->SetAlgSysDataPointer( size, solHelp.GetPointer() );
        singlePDEs_[i]->solVecPrev_  = solVecPrev_;
      }
    }
  }

  // ********************
  //   InitTimeStepping
  // ********************
  void DirectCoupledPDE::InitTimeStepping() {

    PtrParamNode systemNode = FindLinearSystem("direct");

    // Hard Coded
    TS_alg_ = new Newmark( algsys_, systemNode );

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

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->WriteRestart( );
    }
  }

  void DirectCoupledPDE::ReadRestart( UInt &startStep )
  {

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
                                            const Double asteptime ) {

    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->WriteResultsInFile( kstep, asteptime );
    }
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->WriteResultsInFile( kstep, asteptime );
    }

  }


  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void DirectCoupledPDE::InitCoupling(PDECoupling * Coupling)
  {
    isIterCoupled_ = true;
  }
  void DirectCoupledPDE::ResetCoupling()
  {

    iterCoupledCounter_ = 0;
    for (UInt i=0; i<singlePDEs_.GetSize(); i++)
      {
        singlePDEs_[i]->ResetCoupling();
      }
  }

  void DirectCoupledPDE::CalcInputCoupling()
  {
    for (UInt i=0; i<singlePDEs_.GetSize(); i++)
      {
        singlePDEs_[i]->CalcInputCoupling();
      }
  }

  void DirectCoupledPDE::CalcOutputCoupling()
  {
    for (UInt i=0; i<singlePDEs_.GetSize(); i++)
      {
        singlePDEs_[i]->CalcOutputCoupling();
      }
  }




  void DirectCoupledPDE::DefineSolveStep() {

    bool isPiezoHyst  = false;
    bool isMicroPiezo = false;

    // activate direct coupling information
    // and initialize all single pdes
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      // check if single PDE really needs previous solution
      if ( singlePDEs_[i]->BelongsPDE2PiezoHyst() )
        isPiezoHyst = true;

      if (  singlePDEs_[i]->GetName() == "mechanic" ) {
        totalFormulation_ = singlePDEs_[i]->IsTotaFormulation();
      }
      //check for micro-piezo-model
     if ( singlePDEs_[i]->BelongsPDE2MicroPiezo() ) 
        isMicroPiezo = true;
    }

    // check, if we have a mechacou-coupling and have
    // a nonlinear acoustic computation;
    bool acouPDEisNonlinAndCoupled2Mech = false;
    if ( couplings_.GetSize() == 1  &&
         couplings_[0]->GetName() == "acouMechDirect" ) {
         //check if single PDE is acoustic and is nonlinear
         for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
             if ( singlePDEs_[i]->GetName() == "acoustic" &&
                  singlePDEs_[i]->IsNonLin() ) {
                  acouPDEisNonlinAndCoupled2Mech = true;
             }
         }
    }
    if ( acouPDEisNonlinAndCoupled2Mech ) {
      solveStep_ = new SolveStepAcoustic(*this);
    }
    else if ( isPiezoHyst )
      solveStep_ = new SolveStepPiezo(*this);
    else if ( isMicroPiezo )
      solveStep_ = new SolveStepMicroPiezo(*this);
    else
      solveStep_ = new StdSolveStep(*this);
  }

  // *************************
  //   IncorporateReordering
  // *************************
  void DirectCoupledPDE::IncorporateReordering() {


    PdeIdType pdeId   = NO_PDE_ID;
    shared_ptr<EqnMap> eqn;
    StdVector<UInt> newOrder;

    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      pdeId    = singlePDEs_[i]->GetPDEId();
      eqn      = singlePDEs_[i]->GetEqnMap();
      algsys_->GetReordering( pdeId, newOrder );
//       if( newOrder == NULL ) {
//         std::cerr << "performing no reordering!";
//       }

      eqn->ReorderMapping( newOrder );
    }
  }

}
