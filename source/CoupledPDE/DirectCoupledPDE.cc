// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <algorithm>

#include "DirectCoupledPDE.hh"
#include "BasePairCoupling.hh"


// include solveStep drivers
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/SolveSteps/SolveStepHyst.hh"
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
  DirectCoupledPDE::DirectCoupledPDE( Grid *aptgrid, PtrParamNode paramNode,
          PtrParamNode infoNode,
          shared_ptr<SimState> simState,
          Domain* domain)
  
  : StdPDE( aptgrid, paramNode, infoNode, simState, domain ) {
    //		std::cout << "DirectCoupledPDE - Constructor" << std::endl;
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
    //		std::cout << "DirectCoupledPDE - SetSinglePDEs" << std::endl;
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
    olasInfo_ = myInfo_->Get("OLAS")->Get(pdename_);
    
  }
  
  
  // ****************
  //   SetCouplings
  // ****************
  void DirectCoupledPDE::SetCouplings( const StdVector<BasePairCoupling*> &couplings ) {
    couplings_ = couplings;
  }
  
  // ********
  //   Init
  // ********
  void DirectCoupledPDE::Init( UInt sequenceStep )
  {
    //		std::cout << "DirectCoupledPDE - Init" << std::endl;
    sequenceStep_ = sequenceStep;
    
    infoNode_ = myInfo_->Get("PDE")->Get("directCoupledPDE", ParamNode::APPEND);
    
    // Create algebraic system and pass it to SinglePDEs
    isComplex_ = IsComplex();
    algsys_ = new AlgebraicSys(olasNode_, olasInfo_, isComplex_);
    solStrat_ = algsys_->GetSolStrategy();
    
    // ----------------------------
    //  Detection of analysis type
    // ----------------------------
    
    analysistype_ = domain_->GetSingleDriver()->GetAnalysisType( );
    
    // Create new Assemble object
    assemble_ = new Assemble( algsys_, analysistype_, this->mp_, 
            myInfo_->GetRoot() );
    
    
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
      
      // Initialize all SinglePDEs (read domains, materials,
      // define primary results)
      singlePDEs_[i]->Init_Stage1( sequenceStep, infoNode_);
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
      couplings_[i]->Init( sequenceStep_ );
    }
    
    // Perform stage 2 initialization (boundary conditions, integerators)
    for( UInt i = 0; i < singlePDEs_.GetSize(); ++i ) {
      singlePDEs_[i]->Init_Stage2();
    }
    
    // Finalize initialization of SinglePDEs (i.e. FeSpaces, FeFunctions, Time stepping etc.)
    for( UInt i = 0; i < singlePDEs_.GetSize(); ++i ) {
      singlePDEs_[i]->Init_Stage3();
    }
    
    // Initialize all Coupling Objects
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->FinalizeInit( );
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
    
    //     std::cout << "DirectCoupledPDE: pde1_->IsNonLin()? " << singlePDEs_[0]->IsNonLin() << std::endl;
    //     std::cout << "DirectCoupledPDE: pde2_->IsNonLin()? " << singlePDEs_[1]->IsNonLin() << std::endl;
    //     std::cout << "this.IsNonLin()? " << this->nonLin_ << std::endl;
    this->nonLin_ = singlePDEs_[0]->IsNonLin()||singlePDEs_[1]->IsNonLin();
    //     std::cout << "this.IsNonLin()? " << this->nonLin_ << std::endl;
    
    // check for hysteresis; has to be done before solve step gets defined
    isHysteresis_ = false;
		for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
			if((singlePDEs_[i]->IsHysteresis()) && (isHysteresis_ == false) )  {
        //				std::cout << "Hysteresis found!" << std::endl;
				isHysteresis_ = true;
				hysteresisCoefs_ = singlePDEs_[i]->hysteresisCoefs_;
        // we should not need hysteresisPolarization/Magnetization as those 
        // are just the output quantities (i.e. there is nothing to compute from them)
        //				hysteresisPolarization_ = singlePDEs_[i]->hysteresisPolarization_;
        //				hysteresisStrain_ = singlePDEs_[i]->hysteresisStrain_;
			} else if((singlePDEs_[i]->IsHysteresis()) && (isHysteresis_ == true) ) {
				EXCEPTION("Only one of the two pdes should be hysteretic");
			}
    }   
    
    
    // define solveStep-driver
    DefineSolveStep();
    
    // Pass SolveStep object to all single pdes
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      
      singlePDEs_[i]->solveStep_ = solveStep_;
    }
    
    isIterCoupled_ = false;
    
    for ( UInt i = 0; i < singlePDEs_.GetSize(); i++ ) {
      isIterCoupled_ |= singlePDEs_[i]->IsIterCoupled();
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
  
  
  
  
  void DirectCoupledPDE::WriteResultsInFile(const UInt kstep,
          const Double asteptime ) {
    
    for (UInt i=0; i<singlePDEs_.GetSize(); i++) {
      singlePDEs_[i]->WriteResultsInFile( kstep, asteptime );
    }
    for (UInt i=0; i<couplings_.GetSize(); i++) {
      couplings_[i]->WriteResultsInFile( kstep, asteptime );
    }
    
  }
  
  void DirectCoupledPDE::DefineSolveStep() {
    if(isHysteresis_){
			solveStep_ = new SolveStepHyst(*this);
		} else {
			solveStep_ = new StdSolveStep(*this);
		}
  }
	
	
}
