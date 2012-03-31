// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <algorithm>

#include "DirectCoupledPDE.hh"
#include "BasePairCoupling.hh"


// include solveStep drivers
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"

// include PDE classes
#include "PDE/SinglePDE.hh"


#include "Domain/Domain.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TransientDriver.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultHandler.hh"

#include "OLAS/algsys/AlgebraicSys.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  DirectCoupledPDE::DirectCoupledPDE( Grid *aptgrid, PtrParamNode paramNode )

    : StdPDE( aptgrid, paramNode ) {
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

    // Delete BasePairCoupling objects
    for ( UInt i = 0; i < couplings_.GetSize(); i++ ) {
      delete couplings_[i];
    }
    couplings_.Clear();

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

    // create composed name of PDE
    for ( UInt i = 0; i < singlePDEs_.GetSize()-1; i++ ) {
      pdename_ += singlePDEs_[i]->GetName();
      pdename_ += "-";
    }
    pdename_ += singlePDEs_[singlePDEs_.GetSize()-1]->GetName();
    
    // determine unique OLAS id 
    // we collect all ids of the single PDEs and check if they are unique
    std::set<std::string> olasIds;    
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      olasIds.insert( singlePDEs_[i]->olasNode_->Get("id")->As<std::string>() );
    }
    
    if( olasIds.size() > 1 ) {
      EXCEPTION( "Single field PDEs have different algebraic systems assigned." );
    }
    
    olasNode_ = singlePDEs_[0]->olasNode_;
    olasInfo_ = info->Get("OLAS")->Get(pdename_);
    
  }


  // ****************
  //   SetCouplings
  // ****************
  void DirectCoupledPDE::SetCouplings( const StdVector<BasePairCoupling*> &couplings ) {
    couplings_ = couplings;
  }

  // ****************
  //   SetInitial conditions
  //   from the pde members
  // ****************

  void DirectCoupledPDE::SetInitialCondition() {

    REFACTOR;
//    Vector< Double > aux;
//    aux.Init();
//
//    // Construct the initial solution vector
//    shared_ptr<FeSpace> feSpace;
//    UInt singleUnknowns=0;
//    UInt lastIndex=0;
//
//
//    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
//
//      //------------------------------------------------
//      // get the id of the this pde
//      //pdeId = singlePDEs_[i]->GetPDEId();
//      //------------------------------------------------
//
//      // the PDEs members haven't set up their initial conditions yet
//      // -> set the initial conditions of this pde
//      singlePDEs_[i]->SetInitialCondition();
//      if(singlePDEs_[i]->IsSetInitialCondition()==true){
//        this->isSetInitialCondition_=true;
//      }
//
//
//      //------------------------------------------------
//      // get the number of unknowns of this pde
//      singleUnknowns = singlePDEs_[i]->GetNumPdeEquations();
//
//      // check setup of linear system
//      if(singlePDEs_[i]->usePenalty_==false){
//        // uses elimination of Inhomogeneous DBC
//        //std::cout << "Num of Inhomogeneous DBC = "<< eqn->GetNumInHomDirichletEqns () << std::endl;
//        singleUnknowns -= eqn->GetNumInHomDirichletEqns();
//        singleUnknowns -= eqn->GetNumInHomDirichletFileEqns();
//      }
//      //------------------------------------------------
//
//      // init the aux vector
//      // **** it is supoussed that i=pdeID ****
//      // -> the order of the solution vector is the same as ordering of pdeID
//      for (UInt ii = lastIndex; ii < lastIndex+singleUnknowns; ii++) {
//        //aux[ii]=singlePDEs_[i]->getInitialCondition();
//        aux.Push_back(singlePDEs_[i]->getInitialCondition());
//      }
//
//      lastIndex+=singleUnknowns;
//
//    }
//
//    // now we have our solution vector initialized
//    //std::cout << "\n al final aux = "<< aux.Serialize() << std::endl;
//
//    if(this->IsSetInitialCondition()==true){
//
//
//      // save the initial solution vector into algsys
//      algsys_->InitSol(aux);
//
//      // save the initial vector in each pde solution vector
//      SaveSolution(aux.GetPointer(), aux.GetSize());
//
//      // save the initial solution vector into algsys
//      algsys_->InitSol( aux );
//
//    }

  }




  // ********
  //   Init
  // ********
  void DirectCoupledPDE::Init( UInt sequenceStep )
  {
    sequenceStep_ = sequenceStep;

    infoNode_ = info->Get("PDE")->Get("directCoupledPDE", ParamNode::APPEND);
    infoNode_->Get(ParamNode::HEADER)->Get("sequeceStep")->SetValue(sequenceStep);

    
    // Create algebraic system and pass it to SinglePDEs
    algsys_ = new AlgebraicSys(olasNode_, olasInfo_);
    solStrat_ = algsys_->GetSolStrategy();
    
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
    // this problem (e.g. register each pde with FeFctIdType and
    // maximum time derivative order), we should change this here!
    assemble_ = new Assemble( algsys_, analysistype_, 2 );
    

//     // Initialize timestepping
//     if ( analysistype_ == BasePDE::TRANSIENT ) {
//       InitTimeStepping();
//     }

    // activate direct coupling information
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->SetDirectCoupling();
    }
    
    // activate direct coupling information
    // and initialize all single pdes
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->algsys_ = algsys_;
      singlePDEs_[i]->assemble_ = assemble_;
      singlePDEs_[i]->solStrat_ = solStrat_;
      
      // Initialize all SinglePDEs
      singlePDEs_[i]->Init( sequenceStep, infoNode_);
    }

    // Collect all feFunctions defined in single PDEs
    for( UInt i = 0; i < singlePDEs_.GetSize(); ++i ) {
      
     feFunctions_.insert( singlePDEs_[i]->feFunctions_.begin(),
                          singlePDEs_[i]->feFunctions_.end() );
     rhsFeFunctions_.insert(singlePDEs_[i]->rhsFeFunctions_.begin(),
                            singlePDEs_[i]->rhsFeFunctions_.end() );
    }

    // Initialize all Coupling Objects
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->SetAlgSys( algsys_ );
      couplings_[i]->solStrat_ = solStrat_;
      couplings_[i]->SetAssemble( assemble_ );
      couplings_[i]->Init( sequenceStep_, infoNode_ );
    }

    // Finalize initialization of SinglePDEs (i.e. FeSpaces, FeFunctions, Time stepping etc.)
    for( UInt i = 0; i < singlePDEs_.GetSize(); ++i ) {
     singlePDEs_[i]->FinalizeInit();
    }
    
    // Print list of defined integrators of assembly object
    assemble_->ToInfo(infoNode_->Get(ParamNode::HEADER)->Get("integrators"));
    
    // Collect all feFunctions defined in BasePairCouplings
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      feFunctions_.insert( couplings_[i]->feFunctions_.begin(),
                           couplings_[i]->feFunctions_.end() );
      rhsFeFunctions_.insert(couplings_[i]->rhsFeFunctions_.begin(),
                             couplings_[i]->rhsFeFunctions_.end() );
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



  // ********************
  //   InitTimeStepping
  // ********************
  void DirectCoupledPDE::InitTimeStepping() {

    // Pass time time stepping init to the single pdes
    for (UInt i=0, n=singlePDEs_.GetSize(); i<n; i++) {
      singlePDEs_[i]->InitTimeStepping();
    }
  }

  
  void DirectCoupledPDE::UpdateToSolStrategy() {
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->InitTimeStepping();
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
    solveStep_ = new StdSolveStep(*this);
  }


}
