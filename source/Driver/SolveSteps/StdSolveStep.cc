#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include "Utils/mathParser/mathParser.hh"
#include "StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Results/BaseResults.hh"
#include "Utils/EvalIntegrals/BiotSavart.hh"
#include "Utils/Timer.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
//#include "MatVec/SingleVector.hh"
#include "Domain/Results/MHTimeFreqResult.hh"


namespace CoupledField {
  
  // declare logging stream
  DEFINE_LOG(stdsolvestep, "stdsolvestep")
  
  StdSolveStep::StdSolveStep(StdPDE & apde)
  :BaseSolveStep(),
          PDE_(apde)
  {
    
    pdename_      = PDE_.GetName();
    isaxi_        = PDE_.GetIsaxi();
    subdoms_      = PDE_.GetRegions();
    materialData_ = PDE_.GetMaterialData();
    ptgrid_       = PDE_.GetGrid();
    algsys_       = PDE_.GetAlgSys();
    assemble_     = PDE_.GetAssemble();
    solStrat_     = algsys_->GetSolStrategy();
    
    results_      = PDE_.GetResultInfos();
    couplingIter_ = 0;
    solutionLimit_ = NO_SOLUTION_TYPE;
    
    // copy FE functions of PDE
    feFunctions_ = PDE_.GetFeFunctions();
    rhsFeFunctions_ = PDE_.GetRhsFeFunctions();
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator it;
    it = feFunctions_.begin();

    // Copy vectors FE functions in SBM-vector for communication
    // with OLAS and time stepping
    solVec_.SetSize( feFunctions_.size() );
    rhsVec_.SetSize( feFunctions_.size() );

    for( ; it != feFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      FeFctIdType id = ptFct->GetFctId();
      // here the solution vector is filled with pointers from
      // the FE function. Therefore setting the solVec_ in
      // algsys_->GetSolutionVal(solVec_) automatically fills the
      // SingleVector in FE function
      solVec_.SetSubVector(ptFct->GetSingleVector(), id);
    }
    //pos = 0;
    it = rhsFeFunctions_.begin();
    for( ; it != rhsFeFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      FeFctIdType id = ptFct->GetFctId();
      // here the solution vector is filled with pointers from
      // the FE function. Therefore setting the solVec_ in
      // algsys_->GetSolutionVal(solVec_) automatically fills the
      // SingleVector in FE function
      rhsVec_.SetSubVector(ptFct->GetSingleVector(), id);
    }

    // Make sure to have both vectors as "weak" vectors,
    // as the feFunctions themselves are responsible for
    // creation and destruction.
    solVec_.SetOwnership(false);
    rhsVec_.SetOwnership(false);
    
    // nonlinear parameters
    incStopCrit_ = 1e-2;
    residualStopCrit_ = 1e-3;
    nonLinMaxIter_ = 10;
    minValidValue_ = -std::numeric_limits<double>::max();// = DBL_MAX;
    maxValidValue_ = std::numeric_limits<double>::max();// = DBL_MAX;
    
    nonLin_                 = PDE_.IsNonLin();
    nonLinMaterial_         = PDE_.IsNonLinMaterial();
    nonLinTotalFormulation_ = PDE_.IsTotalNonLinFormulation();
    isHyst_                 = PDE_.IsHysteresis();
    regionNonLinTypes_      = PDE_.GetNonLinRegionTypes();
    
    startStep_ = 1;
    
    // set entry type of SBM vector
    //oldRhsLinVal_ = SBM_Vector(BaseMatrix::DOUBLE);
    //tmpOldRhsLinVal_ = SBM_Vector(BaseMatrix::DOUBLE);
    //DeltaRhsLinVal_ =  SBM_Vector(BaseMatrix::DOUBLE);
    RhsLinVal_ = rhsVec_;
    
    //    std::cout << "StdSolveStep - Constructor: " << std::endl;
    //    std::cout << "pdename_: " << pdename_ << std::endl;
    //    std::cout << "isHyst_: " << isHyst_ << std::endl;
        
    // In the end, read nonlinear data from xml-file
    if( nonLin_ || nonLinMaterial_ ) {
      ReadNonLinData();
    }
    
    mHandle_ = PDE_.GetDomain()->GetMathParser()->GetNewHandle();
    mParser_ = PDE_.GetDomain()->GetMathParser();
    mParser_->SetExpr(mHandle_,"step");
  }
  
  
  //! Destructor
  StdSolveStep::~StdSolveStep() {
    logFile_.close();
    mParser_->ReleaseHandle(mHandle_);
  }
  
  
  // ======================================================
  // STATIC SOLVING SECTION
  // ======================================================
  
  void StdSolveStep::PreStepStatic( ) {
    
    // init RHS at this place, because e.g. forces of other PDEs are added
    // to RHS afterwards
    algsys_->InitRHS();
    
  }
  
  
  void StdSolveStep::PostStepStatic() {
    
  }
  
  
  void StdSolveStep::SolveStepStatic() {
    
    if (nonLin_) {
      StepStaticNonLin();
    }
    else {
      StepStaticLin();
    }
  }
  
  
  void StdSolveStep::StepStaticLin() {
    
    assemble_->AssembleMatrices();
    
    // The RHS-sources and boundary conditions
    // have to be reassembled each time
    assemble_->AssembleLinRHS();
    // Set special RHS Values
    PDE_.SetRhsValues();
    
    PDE_.SetBCs();
    
    // store rhs vector back to algsys
    algsys_->GetRHSVal(rhsVec_);
    
    // Only if the matrices have changed (e.g. due to updated lagrangian
    // formulation) the system matrix has to be rebuild
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( NO_FCT_ID, 
              matrix_factor_[NO_FCT_ID] );
    }
    
    // Check if the AMG-framework is used (if so, we have
    // to gather some geometry information at this point)
    if(algsys_->UseAMG() ){
      // only works for elimination
      if( solStrat_->UseDirichletPenalty() ) EXCEPTION("AMG only works for Dirichlet elimination!");
      PDE_.SetGeomInfo();
      algsys_->BuildAMGAuxMatrix();
    }
    
    
    // Incorporate Boundary conditions and
    // recalc the preconditioner eventually
    algsys_->BuildInDirichlet();
    
    if( assemble_->IsMatrixUpdated() ) {
      // get rid of unnecessary zeros (if applicable)
      if ( GetRidOfZerosActive() ) {
        algsys_->GetRidOfZeros(getRidOfZerosTol_);
      }


      algsys_->SetupPrecond();
      
      algsys_->SetupSolver();
    }
    
    // Solve problem
    algsys_->Solve();

    // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
    // afterwards, we notify the solver that the matrix pattern might change again in the next step
    if (GetRidOfZerosActive()) {
      algsys_->RestoreSysMat();
      algsys_->GetSolver()->SetNewMatrixPattern();
    }
    
    // Get the solution and store it
    // Since the entries of solVec_ are pointers to the SingleVector
    // of the FE function, it automatically inserts the values there
    algsys_->GetSolutionVal(solVec_);
  }
  
  
  void StdSolveStep::StepStaticNonLin()
  {
    LOG_DBG(stdsolvestep) << "SSNL: method=" << nonLinMethod_ << " -> " << NonLinMethodTypeEnum.Parse(nonLinMethod_);
    if(NonLinMethodTypeEnum.Parse(nonLinMethod_) == FIXEDPOINT && lineSearch_ == "none")
    {
      StepStaticNonLinFixedPointSimple();
      return;
    }

    // this seems to be a newton implemention, on likely the method setting is ignored.
    // TODO make really clear what is perfomed and change the message accordingly
    nonLinInfo_->Get("implementation")->Get("method")->SetValue("newton");
    if(nonLinMethod_ != "newton")
      nonLinInfo_->SetWarning("selected nonlinear method is '" + nonLinMethod_ + "' but presumable 'newton' is perfomed");

    bool performOneMoreStep;
    bool isNewton = false;

    static_non_lin_step_timer_->Start();
    
    SBM_Vector solInc(BaseMatrix::DOUBLE);
    
    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_; // == current solution
    // =================================
    //  Outer loop: Multilevel strategy 
    // =================================
    UInt numLevels = solStrat_->GetNumSolSteps();

    for( UInt iLevel = 0; iLevel < numLevels; ++iLevel )
    {
      LOG_DBG(stdsolvestep) << "SSNL: iLevel=" << iLevel << " of " << numLevels;
      // update the current solution step in a multilevel approach and
      // inform PDEs (containing the FeSpaces), as well as the AlgebraicSystem
      solStrat_->SetActSolStep(iLevel + 1);
      ReadNonLinData();

      if((lineSearch_ != "none") && (lineSearch_ != "minEnergy")){
        EXCEPTION("The selected line search method is currently only available for energy-based hysteresis");
      }
      
      PDE_.UpdateToSolStrategy();
      algsys_->UpdateToSolStrategy();

      // set the boundary conditions
      PDE_.SetBCs();
      
      //perform the load-steps
      Double loadFactor = 1.0;
      PDE_.GetInfoNode()->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactor);
      
      // setup right hand side
      algsys_->InitRHS();
      Double RhsLinL2Norm = SetLinRHS(loadFactor);
      
      // set iteration counter
      UInt iterationCounter=0;
      
      // =================================
      //  Inner nonlinear loop 
      // =================================
      do {
        iterationCounter++;
        LOG_DBG(stdsolvestep) << "SSNL: ic=" << iterationCounter;
        if ( lineSearch_ != "none" || iterationCounter == 1) {
          //add linear right hand side
          algsys_->InitRHS(RhsLinVal_);
          
          // if the RHS depends on the nonlinearity, we have to re-assemble it
          if( assemble_->IsRhsSolDependent()) {
            assemble_->AssembleNonLinRHS();
          }
          
          // setup the matrices - Be ware that not every PDE actually knows the stiffness matrix!
          // TODO imporove this, e.g. by checking there is actually Newton behind
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);

          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(SYSTEM,solVec_,true);
          solVec_.ScalarMult(-1.0);
        }
        
        
        // assemble Newton bilinear forms
        isNewton = true;
        assemble_->AssembleMatrices(isNewton);
        
        //compute effective matrix
        algsys_->ConstructEffectiveMatrix(NO_FCT_ID, matrix_factor_[NO_FCT_ID]);

        algsys_->BuildInDirichlet();

        // get rid of unnecessary zeros (if applicable)
        if ( GetRidOfZerosActive() ) {
          algsys_->GetRidOfZeros(getRidOfZerosTol_);
        }

        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        
        bool setIDBC = false;
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          setIDBC = true;
        
        algsys_->Solve(setIDBC);

        // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
        // afterwards, we notify the solver that the matrix pattern might change again in the next step
        if (GetRidOfZerosActive()) {
          algsys_->RestoreSysMat();
          algsys_->GetSolver()->SetNewMatrixPattern();
        }

        // new solution is only an increment of the full solution =============
        // Since the entries of solVec_ are pointers to the SingleVector
        // of the FE function, it automatically inserts the values there
        algsys_->GetSolutionVal( solInc, setIDBC );
        //compute norms (residual and incremenal ones)
        Double residualL2Norm = 0.0;
        Double etaLineSearch  = 1.0;
        if ( lineSearch_ == "none" || iterationCounter == 1) {
          //to incooperate the inhomog. Dirichlet BCs we need a full
          //step for the first iteration

          actSol.Add(1.0, solInc);
          // store the new solution
          solVec_ = actSol;
          
          //=================compute residual norm
          algsys_->InitRHS(RhsLinVal_);
          // if the RHS depends on the nonlinearity, we have to re-assemble it
          if( assemble_->IsRhsSolDependent()) {
            assemble_->AssembleNonLinRHS();
          }
          
          // setup the matrices
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);
          
          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(SYSTEM,solVec_,true);
          solVec_.ScalarMult(-1.0);
          
          //get RHS vector
          SBM_Vector actRHS(BaseMatrix::DOUBLE);
          algsys_->GetRHSVal( actRHS );
          
          // calculation of residual error =======================================
          residualL2Norm = actRHS.NormL2();
          //std::cout << "ResAbsolut: " << residualL2Norm << std::endl;
        }
        else {
          // do line search
          residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);
          
          // store the new solution
          solVec_ = actSol;
        }
        
        // calculation relative residual error ====================================
        Double residualErr;
        if ( RhsLinL2Norm > 1.0 )
          residualErr = residualL2Norm / RhsLinL2Norm;
        else
          residualErr = residualL2Norm;
        
        // calculate incremental error ========================================
        Double incrementalErr;
        Double solIncrL2Norm = solInc.NormL2();
        Double actSolL2Norm  = actSol.NormL2();
        
        if ( actSolL2Norm )
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else {
          incrementalErr = solIncrL2Norm;
          WARN("Zero solution vector!! ");
        }

        OutputNonLinIterInfo(pdename_, iLevel+1, iterationCounter, residualErr, incrementalErr, etaLineSearch, PDE_.IsIterCoupled() ? couplingIter_ : -1);

        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep = (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
        
        LOG_DBG(stdsolvestep) << "SSNL: re=" << residualErr << " ie=" << incrementalErr << " poms=" << performOneMoreStep << " maxit=" << nonLinMaxIter_;
        
        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                  << "' in step no '" << iLevel+1 
                  << "' at iteration '" << iterationCounter 
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }
        domain->GetInfoRoot()->ToFile();
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
      
    } // loop over levels
    static_non_lin_step_timer_->Stop();
  }

  void StdSolveStep::StepStaticNonLinFixedPointSimple()
  {
    static_non_lin_step_timer_->Start();

    nonLinInfo_->Get("implementation")->Get("method")->SetValue(NonLinMethodTypeEnum.ToString(FIXEDPOINT));
    nonLinInfo_->Get("implementation")->Get("linesearch")->SetValue("none");
    assert(nonLinMethod_ == NonLinMethodTypeEnum.ToString(FIXEDPOINT));

    SBM_Vector last_sol(solVec_.GetEntryType()); // to compute the diff norm
    last_sol = solVec_; // initiate size and this stuff

    assert(solStrat_->GetNumSolSteps() == 1);

    solStrat_->SetActSolStep(1);
    ReadNonLinData();

    double diff = -1.0;
    bool converged = false;
    for(unsigned int i = 0; i < nonLinMaxIter_ && !converged; i++)
    {
      assemble_->AssembleMatrices();

      algsys_->InitRHS(); // if we don't do this we would sum up over the iterations

      if(assemble_->IsRhsSolDependent())
        assemble_->AssembleNonLinRHS();

      // Set special RHS Values
      PDE_.SetRhsValues();

      PDE_.SetBCs();

      // store rhs vector back to algsys
      algsys_->GetRHSVal(rhsVec_);

      //LOG_DBG2(stdsolvestep) << "SSNLRHS i=" << i << " rhs=" << rhsVec_.ToString();

      // Only if the matrices have changed (e.g. due to updated lagrangian
      // formulation) the system matrix has to be rebuild
      if( assemble_->IsMatrixUpdated() ) {
        algsys_->ConstructEffectiveMatrix( NO_FCT_ID, matrix_factor_[NO_FCT_ID] );
      }

      if(algsys_->UseAMG() ){ // copy and paste from StepStaticLin()
        if(solStrat_->UseDirichletPenalty()) EXCEPTION("AMG only works for Dirichlet elimination!");
        PDE_.SetGeomInfo();
        algsys_->BuildAMGAuxMatrix();
      }

      // Incorporate Boundary conditions and
      // recalc the preconditioner eventually
      algsys_->BuildInDirichlet();

      if(assemble_->IsMatrixUpdated()) {
        // get rid of unnecessary zeros (if applicable)
        if ( GetRidOfZerosActive() ) {
          algsys_->GetRidOfZeros(getRidOfZerosTol_);
        }

        algsys_->SetupPrecond();
        algsys_->SetupSolver();
      }

      algsys_->Solve();

      // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
      // afterwards, we notify the solver that the matrix pattern might change again in the next step
      if (GetRidOfZerosActive()) {
        algsys_->RestoreSysMat();
        algsys_->GetSolver()->SetNewMatrixPattern();
      }

      algsys_->GetSolutionVal(solVec_, true); // true=idbc

      LOG_DBG(stdsolvestep) << "SSNLRHS i=" << " diff=" << solVec_.NormL2(last_sol);
      LOG_DBG2(stdsolvestep) << "SSNLRHS i=" << " actSol=" << last_sol.ToString();
      LOG_DBG2(stdsolvestep) << "SSNLRHS i=" << " sol=" << solVec_.ToString();

      diff = solVec_.NormL2(last_sol) / solVec_.NormL2();

      if(diff <= incStopCrit_)
        converged = true;

      OutputNonLinIterInfo(pdename_, 1, i, -1.0, diff, -1.0);

      last_sol = solVec_;
    }
    static_non_lin_step_timer_->Stop();

    if(!converged && abortOnMaxIter_)
      EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
              << "'within '" << nonLinMaxIter_ << "iterations"
              << "\n ==> incremental error: " << diff);
  }

  // ======================================================
  // Solve Step Transient SECTION
  // ======================================================
  


  void StdSolveStep::InitTimeStepping(){
    //also initialize vectors for the time stepping scheme
    
    stageRHS_.Resize(feFunctions_.size());
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    UInt rhsSize = 0;
    
    //reserve memory for the rhs
    for(fncIt = feFunctions_.begin(); fncIt != feFunctions_.end();++fncIt ){
      rhsSize = fncIt->second->GetSingleVector()->GetSize();
      FeFctIdType id = fncIt->second->GetFctId();
      stageRHS_.SetSubVector(new Vector<Double>(),id);
      stageRHS_.GetPointer(id)->Resize(rhsSize);
    }
    stageRHS_.Init();
  }
  
  void StdSolveStep::PreStepTrans() {
    // Update moving ncInterfaces as needed
    ptgrid_->MoveNcInterfaces();
    
    // due to coupling-pdes, the RHS has to be initialized BEFORE
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();
    PDE_.FinilizeBeforTimeStep();
  }
  
  
  void StdSolveStep::SolveStepTrans() {
    if ( isHyst_ ){
      EXCEPTION("Time stepping for hysteresis no longer implemented in stdsolvestep");
    }
    // do a time step with hysteretic behaviour;
    // TODO ldomenig: isHyst_ is 1/0 to determine if there are any hysteretic regions
    else if (isHyst_){
      StepTransHyst();
    }
    //currently not supported
    //  else if ( nonLin_ && nonLinMaterial_ ) {
    //      StepTransNonLinMaterial();
    //    }
    // do a nonlinear time step
    else if (nonLin_ ){
      if ( nonLinTotalFormulation_ ){
        StepTransNonLinTotal();
      }
      else
        StepTransNonLin();
    }
    // do a linear time step
    else {
      StepTransLin();
    }
  }
  
  
  void StdSolveStep::StepTransLin() {
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    // iterator over solution types and corresponding feFunctions
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    // map to determine the number of temporal derivatives for each matrix types and an iterator for it
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;

    // for iterative coupling and non linear solutions we need to store the initial glmVec
    // for linear transient simulations updatePredictor is always true storeInitialIterGlmVec should be always false
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 ); 
    bool storeInitialIterGlmVec = ( PDE_.IsIterCoupled() == true && couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor,storeInitialIterGlmVec);
    }

    // loop over the number of solver stages of the scheme
    for(UInt i=0;i<numStages;i++){
      rhsVec_.Init();
      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
        FeFctIdType fctId = fncIt->second->GetFctId();
        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),fctId);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);
      algsys_->InitRHS();
      //account for RHS
      assemble_->AssembleLinRHS();
      //Set special RHS Values
      PDE_.SetRhsValues();
      // store rhs vector back to PDE 
      algsys_->GetRHSVal(rhsVec_);

      // assemble matrices...
      // if we want to use static condensation we have to perform the timestepping on element level
      if(algsys_->UseStaticCondensation()){
        matrix_factor_.clear();
        // get timeintegration factors and store them in a map
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()->AddMatFactors(i,matrices,matrix_factor_[fctId]);
        }
        // assemble all matrix parts (needed to update rhs) and also calculate complete system matrix
        assemble_->AssembleMatrices_CondTrans(false,i,matrix_factor_);
        if(assemble_->IsMatrixUpdated()){
          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
            FeFctIdType fctId = fncIt->second->GetFctId();
            // we need to call this function
            // but as we already have summed up all entries on element level we just have to give
            // this function a matrix_factor_ list which contains only 0 else
            std::map<FEMatrixType,Double> zero_factors;
            zero_factors[MASS] = 0.0;
            zero_factors[MASS_UPDATE] = 0.0;
            zero_factors[STIFFNESS] = 0.0;
            zero_factors[STIFFNESS_UPDATE] = 0.0;
            zero_factors[DAMPING] = 0.0;
            zero_factors[DAMPING_UPDATE] = 0.0;
            algsys_->ConstructEffectiveMatrix(fctId, zero_factors); //matrix_factor_[fctId]);
          }
        }
      } else {
        assemble_->AssembleMatrices();
        if(assemble_->IsMatrixUpdated()){
          //if AMG is used, rebuild auxiliary matrix
          auxSet_ = false;
          // set system matrix to zero initially, as ConstructEffectiveMatrix only sums up the contributions
          algsys_->InitMatrix(SYSTEM);
          matrix_factor_.clear();
          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
            FeFctIdType fctId = fncIt->second->GetFctId();
            fncIt->second->GetTimeScheme()->AddMatFactors(i,matrices,matrix_factor_[fctId]);
            algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
          }
        }
      }

      // now compute/update the effective right hand side of the system
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
        if(matIt->second < 0)
          continue;
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(fctId));
        }
        algsys_->UpdateRHS(matIt->first,stageRHS_,assemble_->IsMatrixUpdated());
      }

      // Check if the AMG-framework is used (if so, we have
      // to gather some geometry information at this point)
      // needs only be built once, doesn't change over frequency
      if( (algsys_->UseAMG()) && (auxSet_ == false) ){
        // only works for elimination
        if( solStrat_->UseDirichletPenalty() ) EXCEPTION("AMG only works for Dirichlet elimination!");
        PDE_.SetGeomInfo();
        algsys_->BuildAMGAuxMatrix();
        auxSet_ = true;
      }

      // set boundary conditions
      PDE_.SetBCs();
      algsys_->BuildInDirichlet();

      // prepare the solver and preconditioner for the updated system matrix
      if(assemble_->IsMatrixUpdated() || GetRidOfZerosActive()){

        // get rid of unnecessary zeros (if applicable)
        if ( GetRidOfZerosActive() ) {
          algsys_->GetRidOfZeros(getRidOfZerosTol_);
        }

        algsys_->SetupPrecond();
        algsys_->SetupSolver();
      }

      //write out the current glmVec for debugging purposes
      if (IS_LOG_ENABLED(stdsolvestep, dbg3)){
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end(); ++fncIt){
          LOG_DBG(stdsolvestep) <<"PDE name: " << fncIt->second->GetPDE()->GetName() << std::endl;
          LOG_DBG(stdsolvestep) <<"feFctId: " << fncIt->second->GetFctId() << std::endl;
          LOG_DBG(stdsolvestep) <<"actStep: " << this->actStep_ << std::endl;
          LOG_DBG(stdsolvestep) <<"couplingIter: " << this->couplingIter_ << std::endl;
          LOG_DBG(stdsolvestep) <<"GLM vector: " << std::endl;

          for(UInt i=0;i<fncIt->second->GetTimeScheme()->GetSizeGLMVector();i++){
            SingleVector* glmVector = fncIt->second->GetTimeScheme()->GetInitialIterGLMVector(i);
            LOG_DBG(stdsolvestep) << "Index " << i << std::endl;
            LOG_DBG(stdsolvestep) << glmVector->ToString(TS_NONZEROS,"\n") << std::endl;
            LOG_DBG(stdsolvestep) << "Finish GLM Vector" << std::endl;
          }
          LOG_DBG(stdsolvestep) << std::endl;

          LOG_DBG(stdsolvestep) << "This is the initial GLM Vector of this time step" << std::endl;
          for(UInt i=0;i<fncIt->second->GetTimeScheme()->GetSizeGLMVector();i++){
            SingleVector* initialIterGlmVector = fncIt->second->GetTimeScheme()->GetInitialIterGLMVector(i);
            LOG_DBG(stdsolvestep) << "Index " << i << std::endl;
            LOG_DBG(stdsolvestep) << initialIterGlmVector->ToString(TS_NONZEROS,"\n") << std::endl;
            LOG_DBG(stdsolvestep) << "Finish initial GLM Vector" << std::endl;
          }
          LOG_DBG(stdsolvestep) << std::endl;
        }
      }
      // solve the system of equations using the defined solver
      algsys_->Solve();

      // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
      // afterwards, we notify the solver that the matrix pattern might change again in the next step
      if (GetRidOfZerosActive()) {
        algsys_->RestoreSysMat();
        algsys_->GetSolver()->SetNewMatrixPattern();
      }
      
     // Since the entries of solVec_ are pointers to the SingleVector
     // of the FE function, it automatically inserts the values there
      algsys_->GetSolutionVal(stageSol);
    }
    //update stage
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end(); ++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }
  }
  
  
  void StdSolveStep::StepTransNonLin() {
    /*!
     * Comments added to better understand what's going on. If you find any errors, please correct.
     */
    /*! Solve step for transient, nonlinear simulations using material non-linearities
     *  Works for std fixpoint iteration as well as for Newton iteration.
     *  Instead of computing the total value u_n+1^k+1 at time step n+1 and iteration k+1,
     *  only the difference deltaU^k+1 = u_n+1^k+1 - u_n+1^k is computed.
     *
     *  This leads to the following iterative scheme (only fixpoint case described here)
     *    n = time step counter; k = iteration counter
     *    K(u) = nonlinear system matrix; f = rhs load
     *    Known from last iteration: u_n+1^k
     *    Aim: compute u_n+1^k+1 such that
     *      a) u_n+1^k+1 - u_n+1^k = deltaU^k+1 < tol        -> small incremental error
     *      b) K( u_n+1^k+1 )*u_n+1^k+1 - f     < tol        -> small residual error
     *
     *    Approach:
     *    1. Build up system:
     *      K( u_n+1^k )*u_n+1^k+1 = f
     *    > K( u_n+1^k )*(u_n+1^k + deltaU^k+1) = f
     *    > K( u_n+1^k )*deltaU^k+1 = f - K( u_n+1^k )*u_n+1^k
     *
     *    2. solve for solution increment deltaU^k+1
     *
     *    3. calculate new solution by updating the old one
     *      u_n+1^k+1 = u_n+1^k + deltaU^k+1
     *
     *    4. check for incremental and residual error
     *
     *    5. continue if one of the errors is too large and the maximal number of iterations
     *        has not been reached yet
     *
     *    Three possible implementations (I and II work well at least for the NLtrans magnetic example from
     *      testsuite, III was not tested in that context)
     *      I. u_n+1^0 = u_n
     *         u_n+1^1 = deltaU^1
     *         u_n+1^2 = deltaU^1 + deltaU^2
     *
     *         This is how it is implemented below. During the first iteration, we calculate
     *         deltaU^1 between u_n and u_n+1^1, but we do NOT add deltaU^1 to u_n
     *         to approximate u_n+1^1.
     *         Instead we set u_n+1^1 = deltaU^1.
     *         This leads to (completely) different evaluation points for the material curve and seems
     *         a little bit odd at first. However, it converges to the same solution as II.
     *
     *      II. u_n+1^0 = u_n, but rhs will NOT be updated during first iteration
     *         -> we solve not for deltaU^1 but for the whole vector u_n+1^1 during the first iteration
     *         -> i.e. we perform one "total" step (like it is done in StepTransNonLinTotal())
     *
     *         K( u_n+1^0 )*u_n+1^1 = f
     *         -> u_n+1^1 will thus contain u_n
     *
     *         u_n+1^2 = u_n + deltaU^1 + deltaU^2
     *         -> all further iterations are standard delta-steps
     *         -> compared to I. we do not evaluate the material tensor only with the summed up increments,
     *            but with the increments added to the previous solution
     *         -> tests showed, that I. and II. will converge to the same result, although the calculated steps
     *            to come there are different
     *
     *      III. u_n+1^0 = u_n
     *         u_n+1^1 = u_n + deltaU^1
     *         u_n+1^2 = u_n + deltaU^1 + deltaU^2
     *         ...
     *         -> solution from old time step is the basis for the new solution and will be kept
     *         -> seems to me the most natural way but it has the problem when the system is solved
     *            for the first time during each time step. To incorporate changing inhomogeneous Dirichlet
     *            Boundary Conditions (IDBC) into the result vector, we would have to replace the IDBC values
     *            from the last time step with the ones from the new time step. A simple addition would lead
     *            to wrong results here
     *         -> Added function ClearIDBCFromSolutionVal to set all IDBC nodes to 0
     *            -> call this function with the OLD solution before solving the system
     *            -> the old solution will have 0 at IDBC nodes so that the new solution can simply be added
     *
     *
     *
     *
     */
    
    mParser_->SetExpr(MathParser::GLOB_HANDLER,"iterationCounter");
    mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", 0);

    LOG_TRACE(stdsolvestep) << "StdSolveStep::StepTransNonLin";
    bool performOneMoreStep;
    bool isNewton = false;
    
    SBM_Vector solInc(BaseMatrix::DOUBLE);
    
    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;
    
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;
    
    UInt pos = 0;
    
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 ); 
    bool storeInitialIterGlmVec = ( couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor,storeInitialIterGlmVec);
    }
    
    for(UInt i=0;i<numStages;i++){
      //do initialization 
      rhsVec_.Init();
      LOG_DBG(stdsolvestep) << "StepTransNonLin: Stage: " << i ;
      
      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        FeFctIdType fncId = fncIt->second->GetFctId();

        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),fncId);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);
      
      
      //initialize solution vector for each stage
      if ( i > 0 )
        actSol = stageSol;
      else{
        //special case of incremental non-linearity, we set the stage vector to the solution vector
        stageSol = actSol;
      }
      
      solVec_  = actSol;
      
      // setup right hand side
      Double loadFactor = 1.0;
      Double RhsLinL2Norm = SetLinRHS(loadFactor);
      
      // set iteration counter
      UInt iterationCounter=0;
      
      do {
        iterationCounter++;
        mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", iterationCounter);

        if ( lineSearch_ != "none" || iterationCounter == 1) {
          //add linear right hand side
          algsys_->InitRHS(RhsLinVal_);
          
          // if the RHS depends on the nonlinearity, we have to re-assemble it
          if( assemble_->IsRhsSolDependent() ) {
            assemble_->AssembleNonLinRHS();
          }
          
          // setup the matrices
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);
          
          //now update RHS according to time stepping
          for(matIt = matrices.begin();matIt != matrices.end();matIt++){
            
            if(matIt->second < 0)
              continue;
            for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
              FeFctIdType fncId = fncIt->second->GetFctId();
              fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(fncId));
            }
            algsys_->UpdateRHS(matIt->first,stageRHS_,true);
          }
          
          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(STIFFNESS,solVec_,true); // we also or only need the updated version
          algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
          solVec_.ScalarMult(-1.0);
        }
        
        //now assemble the Newton bilinear forms
        isNewton = true;
        assemble_->AssembleMatrices(isNewton);
        
        matrix_factor_.clear();
        
        // set system matrix to zero initially, as ConstructEffectiveMatrix only
        // sums up the contributions
        algsys_->InitMatrix(SYSTEM);
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()
          ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
          algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
        }
        // setup the matrices to compute correct error norms
        
        PDE_.SetBCs();
        algsys_->BuildInDirichlet();

        // get rid of unnecessary zeros (if applicable)
        if ( GetRidOfZerosActive() ) {
          algsys_->GetRidOfZeros(getRidOfZerosTol_);
        }

        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        
        // just set inh. Dirichlet BCs for the first iteration
        bool setIDBC = false;
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          setIDBC = true;
        
        algsys_->Solve(setIDBC);

        // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
        // afterwards, we notify the solver that the matrix pattern might change again in the next step
        if (GetRidOfZerosActive()) {
          algsys_->RestoreSysMat();
          algsys_->GetSolver()->SetNewMatrixPattern();
        }

        // if setIDBC is true, solInc will contain the inhom. Dirichlet values
        // Since the entries of solVec_ are pointers to the SingleVector
        // of the FE function, it automatically inserts the values there
        algsys_->GetSolutionVal(solInc, setIDBC );
        
        Double residualL2Norm = 0.0;
        Double etaLineSearch  = 1.0;
        
        //necessary due to inh. Dirichlet BCs!!
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          stageSol.Init();
        
        if ( lineSearch_ == "none" || iterationCounter == 1) {
          //to incooperate the inhomog. Dirichlet BCs we need a full
          //step for the first iteration
          
          /*
           * note: stageSol was initialized to 0.0, so
           *      during the iterations, stageSol will be the sum of all solultion increments but will NOT
           *      contain the old solution itself
           */
          stageSol.Add(1.0, solInc);
          
          solVec_  = stageSol;
          
          //=================compute residual norm
          algsys_->InitRHS(RhsLinVal_);
          // if the RHS depends on the nonlinearity, we have to re-assemble it
          if( assemble_->IsRhsSolDependent()) {
            assemble_->AssembleNonLinRHS();
          }
          
          // setup the matrices with new solution
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);
          
          //now update RHS according to time stepping
          for(matIt = matrices.begin();matIt != matrices.end();matIt++){
            //std::cout << "Matrix: " << matIt->first << std::endl;
            if(matIt->second < 0)
              continue;
            for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
              FeFctIdType fncId = fncIt->second->GetFctId();
              fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(fncId));
            }
            algsys_->UpdateRHS(matIt->first,stageRHS_,true);
          }
          
          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(STIFFNESS,solVec_,true);
          algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
          solVec_.ScalarMult(-1.0);
          
          //get RHS vector
          SBM_Vector actRHS(BaseMatrix::DOUBLE);
          algsys_->GetRHSVal( actRHS );
          
          // calculation of residual error =======================================
          residualL2Norm = actRHS.NormL2();
        }
        else {
          //solVec_  = stageSol;
          residualL2Norm = LineSearch(solInc, stageSol, etaLineSearch,true);
          solVec_  = stageSol;
        }
        
        // calculation of residual error =======================================
        Double residualErr;
        if ( RhsLinL2Norm > 1.0 )
          residualErr = residualL2Norm / RhsLinL2Norm;
        else
          residualErr = residualL2Norm;
        
        // calculate incremental error ========================================
        Double incrementalErr;
        Double solIncrL2Norm = solInc.NormL2();
        Double actSolL2Norm  = stageSol.NormL2();
        
        if ( actSolL2Norm )
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else {
          incrementalErr = solIncrL2Norm;
          //WARN("Zero solution vector!! ");
        }

        OutputNonLinIterInfo(pdename_, PDE_.GetSolveStep()->GetActStep(),iterationCounter, residualErr, incrementalErr, etaLineSearch, PDE_.IsIterCoupled() ? couplingIter_ : -1);

        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep = (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
        
        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                  << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
                  << "' at iteration '" << iterationCounter 
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
      
    } //stages
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator limitFeFctIt;
    limitFeFctIt = feFunctions_.find(solutionLimit_);
    if (limitFeFctIt != feFunctions_.end() ) {
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        FeFctIdType fncId = fncIt->second->GetFctId();
        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
          //const SingleVector * subv = solVec_.GetPointer(pos);
          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(fncId)));
          for (UInt j=0; j < dsubVec.GetSize(); j++) {
            if (dsubVec[j] >= maxValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ << 
                      "' is larger than the allowed maximum limit set in the XML: "
                      << maxValidValue_); 
            }
            if (dsubVec[j] <= minValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ << 
                      "' is smaller than the allowed minimum limit set in the XML: "
                      << minValidValue_); 
            }
          }
        }
      }
    }
    
    //update stage
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      /*
       * here we finally compute the new solution vector
       *  solution_new = solution_old + stage_solutions
       */
      fncIt->second->GetTimeScheme()->FinishStep();
    }
  }
  
  
  void StdSolveStep::StepTransNonLinTotal() {
    
    bool performOneMoreStep;
    bool isNewton;
    Double incrementalErr;
    
    SBM_Vector solNew(BaseMatrix::DOUBLE);
    SBM_Vector diffSol(BaseMatrix::DOUBLE);
    
    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;
    
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;
    
    UInt pos = 0;
    
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
    bool storeInitialIterGlmVec = ( couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor,storeInitialIterGlmVec);
    }
    
    for(UInt i=0;i<numStages;i++){
      //do initialization
      rhsVec_.Init();
      
      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),pos);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);
      
      //initialize solution vector for each stage
      if ( i > 0 )
        actSol = stageSol;
      else{
        //special case of incremental non-linearity, we set the stage vector to the solution vector
        stageSol = actSol;
      }
      
      solVec_  = actSol;
      
      // setup right hand side
      Double loadFactor = 1.0;
      incrementalErr = SetLinRHS(loadFactor);
      
      // set iteration counter
      UInt iterationCounter=0;
      
      do {
        iterationCounter++;
        
        // do matrices: Newton is not working for total formulation!!
        isNewton = false;
        assemble_->AssembleMatrices(isNewton);
        
        // set RHS
        algsys_->InitRHS(RhsLinVal_);
        
        /*
         * set nonlinrhs via assmble and not via setlinrhs(...,true)
         */
        assemble_->AssembleNonLinRHS();
        
        //now update RHS according to time stepping
        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
          if(matIt->second < 0)
            continue;
          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
            fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
          }
          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
        }
        
        // set system matrix to zero initially, as ConstructEffectiveMatrix only
        // sums up the contributions
        matrix_factor_.clear();
        algsys_->InitMatrix(SYSTEM);
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()
          ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
          algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
        }
        
        PDE_.SetBCs();
        algsys_->BuildInDirichlet();

        // get rid of unnecessary zeros (if applicable)
        if ( GetRidOfZerosActive() ) {
          algsys_->GetRidOfZeros(getRidOfZerosTol_);
        }

        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        
        //always set inhomog. Dirichlet BCs
        bool setIDBC = true;
        
        algsys_->Solve(setIDBC); 

        // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
        // afterwards, we notify the solver that the matrix pattern might change again in the next step
        if (GetRidOfZerosActive()) {
          algsys_->RestoreSysMat();
          algsys_->GetSolver()->SetNewMatrixPattern();
        }

        // Since the entries of solVec_ are pointers to the SingleVector
        // of the FE function, it automatically inserts the values there
        algsys_->GetSolutionVal(solNew, setIDBC );
        
        // calculate incremental error ========================================
        diffSol = solNew;
        diffSol.Add( -1.0, actSol);
        Double solIncrL2Norm = diffSol.NormL2();
        Double solNewL2Norm = solNew.NormL2();
        
        if (solNewL2Norm > 1)
          incrementalErr = solIncrL2Norm / solNewL2Norm;
        else
          incrementalErr = solIncrL2Norm;
        
        //just dummy things
        Double etaLineSearch = 1.0;
        Double residualErr = incrementalErr;
        
        OutputNonLinIterInfo(pdename_, 1,iterationCounter, residualErr, incrementalErr, etaLineSearch);

        stageSol = solNew;
        solVec_  = stageSol;
        
        //store new solution
        actSol = solNew;
        
        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
                (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
        
        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                  << "' at iteration '" << iterationCounter 
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }
        
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
      
    } //stages
    
    //update stage
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }
  }
/** direct quasi-Newton formlation
 * Solves for one transient time step by using the direct quasi-Newton formulation.
 * In this formulation the PDE is directly linearized by a Taylor Series expansion.
 * So we come up with a linear algebraic system of the form
 * dF/du * delta_u = -F
 * where F is the residual of the PDE and dF/du is the Jacobian of the PDE.
 * 
 * ### MAGNETOSTATIC PHI-FORMULATION ###
 * This would lead to the following weak form
 * (-dB/dH gradDelta_Phi, gradPhi')_Omega = -(B_prev, gradPhi')_Omega .
 * 
 * ### MAGNETOSTATIC A-FORMULATION ###
 * This would lead to the following weak form
 * (-dH/dB curlDelta_A, curlA')_Omega = -(H_prev, curlA')_Omega - (J_s, A')_Omega .
 * 
 * where dB/dH or dH/dB is the material tensor mu(H) or nu(B) respectively.
 * This system of equations has to be solved until the error norm is small enough.
 * 
 * =================================================================================
 * ### Start: ALGORITHM
 * =================================================================================
 * u = linearFEM() // initial guess with current boundary conditions
 * while(true) // nonlinear iterations
 *   getMaterialTensor() // evaluate hysteresis operator + DFP for mu(H)
 *   solveSystem() // solve linearized system from above
 *   performNewtonStep() // maybe also use a linesearch method
 *   checkStoppingCriteria() // check if error small enough
 *   saveState() // save all needed quantities for next iteration (if necessary)
 * end while
 * 
 * saveState() // save all needed quantities for next time step
 * =================================================================================
 * ### End: ALGORITHM
 * =================================================================================
**/
  void StdSolveStep::StepTransHyst(){
    mParser_->SetExpr(MathParser::GLOB_HANDLER,"iterationCounter");
    mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", 0);

    LOG_TRACE(stdsolvestep) << "StdSolveStep::StepTransHyst";
    bool performOneMoreStep;


    // TODO ldomenig: get initial guess from linear case with current boundary conditions
    // and  previous mu


    // =================================================================================
    // ### nonlinear iterations
    // =================================================================================
    // set iteration counter
    UInt iterationCounter = 0;

    do { // nonlinear iteration start
      iterationCounter++;
      mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", iterationCounter);
      SBM_Vector deltaPhi(BaseMatrix::DOUBLE);




      Double deltaPhiL2Norm = 0.0;
      // calculation of ||Delta_phi||_L2 ==============================================
      deltaPhiL2Norm = deltaPhi.NormL2();


      // check if one more step is necessary ==========================================
      performOneMoreStep = (deltaPhiL2Norm > incStopCrit_);
      if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_){
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                    << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
                    << "' at iteration '" << iterationCounter 
                    << "'.\n ==> delta phi L2norm: " << deltaPhiL2Norm);
      }
      if (performOneMoreStep){
        // TODO ldomenig: save current state for next iteration
        
      }
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_); // nonlinear iteration end

      // TODO ldomenig: save current state for next time step
  } 

  
  // ======================================================
  // Solve Step Harmonic  SECTION
  // ======================================================w
  
  void StdSolveStep::PreStepHarmonic() {
    algsys_->InitRHS();
  }

  void StdSolveStep::SolveStepHarmonic() {
    if ( nonLin_ || solStrat_->IsMultHarm() ) {
      StepHarmonicNonLin();
    }
    else {
      StepHarmonicLin();
    }
  }
  
  
  void StdSolveStep::StepHarmonicLin() {
    //Set special RHS Values
    //std::cout << "Do Apply Loads" << std::endl;
    PDE_.SetRhsValues();
    
    //this has to be done each frequency!
    assemble_->AssembleLinRHS();

    assemble_->AssembleMatrices( );
    PDE_.SetBCs();
    
    // store rhs vector back to PDE
    algsys_->GetRHSVal( rhsVec_ );

    // Where should we get the matrix factors from in a harmonic case?
    // In my opinion this method
    //if( assemble_->IsMatrixUpdated() ) {
    std::map<FEMatrixType,Double> empty;
    algsys_->ConstructEffectiveMatrix(NO_FCT_ID,  empty );
    

    // Check if the AMG-framework is used (if so, we have
    // to gather some geometry information at this point)
    // needs only be built once, doesn't change over frequency
    if( ((algsys_->UseAMG()) && (auxSet_ == false)) ||
            ((algsys_->UseAMG()) && (algsys_->GetAMGType() != AMGType::EDGE)) ){
      // only works for elimination
      if( solStrat_->UseDirichletPenalty() ) EXCEPTION("AMG only works for Dirichlet elimination!");
      PDE_.SetGeomInfo();
      algsys_->BuildAMGAuxMatrix();
      auxSet_ = true;
    }
    
    // Incorporate Boundary conditions and
    // recalc the preconditioner eventually
    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      // get rid of unnecessary zeros (if applicable)
      if ( GetRidOfZerosActive() ) {
        algsys_->GetRidOfZeros(getRidOfZerosTol_);
      }

      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }
    
    algsys_->Solve();

    // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
    // afterwards, we notify the solver that the matrix pattern might change again in the next step
    if (GetRidOfZerosActive()) {
      algsys_->RestoreSysMat();
      algsys_->GetSolver()->SetNewMatrixPattern();
    }

    // Since the entries of solVec_ are pointers to the SingleVector
    // of the FE function, it automatically inserts the values there
    algsys_->GetSolutionVal(solVec_);
    
    if ( adjointSource_ ) {
      //check if adjoint PDE has been solved in case of source localization
      //if yes, we have to multiply the solution with a standard mass matrix
      std::cout << "DO multiply with MASS-matrix" << std::endl;
      //solVec_.Export("sol1.dat",BaseMatrix::MATRIX_MARKET);
      algsys_->InitRHS();
      algsys_->UpdateRHS(AUXILIARY,solVec_,true);
      algsys_->GetRHSVal( solVec_ );
      //solVec_.Export("sol2.dat",BaseMatrix::MATRIX_MARKET);
      //std::cout << "SOL after: \n " << solVec_ << std::endl;
      adjointSource_ = false;
    }
  }

  void StdSolveStep::StepHarmonicNonLin() {

    LOG_TRACE(stdsolvestep) << "LineSearch used: " << lineSearch_ << std::endl;
    // Set some variables
    UInt N = solStrat_->GetNumHarmN();
    UInt M = solStrat_->GetNumHarmM();
    Double bF = solStrat_->GetBaseFreq();
    UInt numFFT = solStrat_->GetNumFFT();
    if(numFFT % 2 != 0){
      EXCEPTION("Please provide a numFFT xml attribute, which is even!");
    }

    //setting the iterationCounter so one can easy check which iteration it is
    mParser_->SetExpr(MathParser::GLOB_HANDLER,"iterationCounter");
    mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", 0);

    bool performOneMoreStep = true;
    // =================================================================================
    //  1) Solve the initial multiharmonic ''linear'' system
    // =================================================================================

    // Perform the load-steps
    Double loadFactor = 1.0;
    PDE_.GetInfoNode()->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactor);

    // setup right hand side
    algsys_->InitRHS();
    // first boolean is flag if nonlinear (first iteration is linear)
    // second boolean is if it's multiharmonic...which is is
    Double RhsLinL2Norm = SetLinRHS(loadFactor, false, true);

    if (IS_LOG_ENABLED(stdsolvestep, dbg3)) std::cout<<"Right Hand Side Linear "<<RhsLinVal_.ToString()<<std::endl;

    // Usually the RhsLinVal_ gets set in the constructor but
    // not in the multiharmonic case. Therefore we set it here.
    // Already done by SetLinRHS()
    //algsys_->GetFullMultiHarmRHSVal(RhsLinVal_);

    // Loop over every frequency and assemble the correct SBM blocks
    AssembleMH(N, M, true);
    // Sets flag that matrix was already assembled. The method CheckNonLinearities
    // redoes this
    assemble_->PostAssemble();


    // Calls method ApplyBC and ApplyLoads in FeFunction
    PDE_.SetBCs();

    // Computation of effective matrix:
    /* NOTE: this is commented because we also include the MASS matrix
             in the SYSTEM matrix, as defined in Assemble::CreateMatrixMap().
             Sometimes having an extra MASS matrix is benefitial, e.g. for exporting and
             comparing different matrix parts, that's why it's still here

       NOTE2: uncommented because we have to set IDBCs
       (Afterwards a different mechanism handles the IDBC)
    */
    std::map<FEMatrixType,Double> empty;
    algsys_->ConstructEffectiveMatrix(NO_FCT_ID,  empty, true );


    // Incorporate Boundary conditions and
    // recalculate the preconditioner eventually
    algsys_->BuildInDirichlet();

    // get rid of unnecessary zeros (if applicable)
    if ( GetRidOfZerosActive() ) {
      algsys_->GetRidOfZeros(getRidOfZerosTol_);
    }

    algsys_->SetupPrecond();
    algsys_->SetupSolver();

    // Solve the linear system
    algsys_->Solve();

    // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
    // afterwards, we notify the solver that the matrix pattern might change again in the next step
    if (GetRidOfZerosActive()) {
      algsys_->RestoreSysMat();
      algsys_->GetSolver()->SetNewMatrixPattern();
    }

    // Get the solution of the initial (linear) multiharmonic system.
    solVecMH_.ResetEntryType(BaseMatrix::EntryType::COMPLEX);
    algsys_->GetFullMultiHarmSolutionVal( solVecMH_, true); //was false, but now we have IDBCs

    if (IS_LOG_ENABLED(stdsolvestep, dbg3)) std::cout<<"SOLUTION OF LINEAR SYSTEM"<<solVecMH_.ToString()<<std::endl;


    // Get actual solution. Usually it is done via actSol = solVec_;
    // but we need the full multiharmonic solution vector
    SBM_Vector actSol(BaseMatrix::COMPLEX);
    actSol = solVecMH_;

    // Create multiharmonic time-frequency object and provide basic information
    MHTimeFreqResult ftRes(N, M, bF, numFFT, PDE_.GetDomain());

    // Evaluate the nonlinearity (transform solution in time domain =>
    // evaluate the curl for the B-field => evaluate BH curve =>
    // transform the nu(t) back to frequency domain nu(harmonic)
    this->EvaluateNonlinearity(ftRes, actSol);

    if (IS_LOG_ENABLED(stdsolvestep, dbg3)){
      std::cout<<"actSol = "<<actSol.ToString()<<std::endl;
    }

    // =================================================================================
    //  2) Solve the full multiharmonic nonlinear system
    // =================================================================================

    // Create new timer object and put it to related info element
    shared_ptr<Timer> timer(new Timer());
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    // We don't have a TwoLevel Strategy!!!
    iter->GetByVal("solStep","value",1,ParamNode::INSERT)->Get("timer")->SetValue(timer);
    timer->Start();

    // As long as we don't have a TwoLevel solution strategy this method
    // does not do anything...
    solStrat_->SetActSolStep(1);

    // This method just reads the nonlinear xml-node
    ReadNonLinData();

    // UpdateToSolStrategy in FeSpaces is actually only used for the TwoLevel
    // solution strategy, where the first level only contains lowest order
    // Hcurl basis functions and the other level higher order basis.
    // In our case this method does not do anything
    PDE_.UpdateToSolStrategy();

    // set iteration counter
    UInt iterationCounter = 0;


    // This will be the incremental solution (deflect vector),
    // meaning \Delta u^{k+1} = u^{k+1} - u^k
    SBM_Vector solInc(BaseMatrix::COMPLEX);

    // ===============================================================
    //  2.1) Nonlinear loop
    // ===============================================================
    lastError_=10000;
    do {
      iterationCounter++;
      mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", iterationCounter);

      // if the RHS depends on the nonlinearity, we have to re-assemble it
      if( assemble_->IsRhsSolDependent() ) {
        EXCEPTION("StdSolveStep::StepHarmonicNonLin() cannot handle solution-dependent"
            "RHSs yet!")
      }

      AssembleMH(N, M);
      // Sets flag that matrix was already assembled. The method CheckNonLinearities re-does this
      assemble_->PostAssemble();

      // Computation of effective matrix:
      /* NOTE: this is commented because we also include the MASS matrix
               in the SYSTEM matrix, as defined in Assemble::CreateMatrixMap().
               Sometimes having an extra MASS matrix is benefitial, e.g. for exporting and
               comparing different matrix parts, that's why it's still here
      */
      //std::map<FEMatrixType,Double> empty;
      //algsys_->ConstructEffectiveMatrix(NO_FCT_ID,  empty, true );


      // set RHS: linear part
      algsys_->InitRHS(RhsLinVal_ );


      if (IS_LOG_ENABLED(stdsolvestep, dbg3)) {
        SBM_Vector resid(BaseMatrix::COMPLEX);
        algsys_->GetFullMultiHarmRHSVal(resid);
        std::cout<<"RHSLINVAL OF STEP "<<iterationCounter<<" = \n"<<RhsLinVal_.ToString()<<std::endl;
        std::cout<<"RESIDUAL VECTOR OF STEP "<<iterationCounter<<" = \n"<<resid.ToString()<<std::endl;
        std::cout<<"MHSOLVEC VECTOR OF STEP "<<iterationCounter<<" = \n"<<solVecMH_.ToString()<<std::endl;
      }

      // This is done because we want to solve the deflect-system:
      // K(u^k) \cdot \Delta u^{k+1} = f - K(u^k) \cdot u^k
      // where f - K(u^k) \cdot u^k gets set as the new rhs in algsys_
      solVecMH_.ScalarMult(-1.0);
      // the boolean has no effect...
      algsys_->UpdateRHS_MultHarm(SYSTEM,solVecMH_,true);
      solVecMH_.ScalarMult(-1.0);

      if (IS_LOG_ENABLED(stdsolvestep, dbg3)) {
        SBM_Vector resid(BaseMatrix::COMPLEX);
        algsys_->GetFullMultiHarmRHSVal(resid);
        std::cout<<"RESIDUAL VECTOR OF STEP "<<iterationCounter<<" = \n"<<resid.ToString()<<std::endl;
        std::cout<<"MHSOLVEC VECTOR OF STEP "<<iterationCounter<<" = \n"<<solVecMH_.ToString()<<std::endl;
      }

      // Incorporate Boundary conditions and recalc the preconditioner and solver
      algsys_->BuildInDirichlet();

      // get rid of unnecessary zeros (if applicable)
      if ( GetRidOfZerosActive() ) {
        algsys_->GetRidOfZeros(getRidOfZerosTol_);
      }

      algsys_->SetupPrecond();
      algsys_->SetupSolver();

      // Solve the deflect system K(u^k) \cdot \Delta u^{k+1} = f - K(u^k) \cdot u^k
      // for the deflect-vector \Delta u^{k+1}
      // DO WE NEED TO CALL IT WITH SETIDBC? No only in the linear iteration.
      // Now its handled in updateRHS_MultHarm()
      algsys_->Solve(false);

      // we store the old (non-optimized) matrix back IMMIDEATELY so that all matrix update operations work again
      // afterwards, we notify the solver that the matrix pattern might change again in the next step
      if (GetRidOfZerosActive()) {
        algsys_->RestoreSysMat();
        algsys_->GetSolver()->SetNewMatrixPattern();
      }

      // Get the incremental solution (deflect vector), second argument is setIDBC
      algsys_->GetFullMultiHarmSolutionVal( solInc, false);


      if (IS_LOG_ENABLED(stdsolvestep, dbg3)) std::cout<<"SOLUTION INCREMENT AT STEP "<<iterationCounter<<" = \n"<<solInc.ToString()<<std::endl;

      // Initialize norms (residual and incremental ones)
      Double residualL2Norm = 0.0;
      Double etaLineSearch  = 1.0;


      // Perform line search to get the 'optimal' eta, which minimizes the residual-norm
      // Meaning: u^{k+1} = u^k + eta * \Delta u^{k+1} in order to minimize
      // the residual r^{k+1} = f - K(u^{k+1}) \cdot u^{k+1} is minimized
      residualL2Norm = LineSearchMultHarm(solInc, actSol, etaLineSearch, ftRes);

      this->EvaluateNonlinearity(ftRes, actSol);

      if (IS_LOG_ENABLED(stdsolvestep, dbg3)) std::cout<<"SOLUTION VECTOR AT STEP "<<iterationCounter<<" = \n"<<actSol.ToString()<<std::endl;

      // Store the new solution u^{k+1}
      // Usually actSol is stored in solVec_ but this is not our full multiharmonic
      // solution vector, therefore we store it in the temporary multiharmonic
      // solution vector solVecMH
      solVecMH_ = actSol;

      if (IS_LOG_ENABLED(stdsolvestep, dbg3)) std::cout<<"SOLUTION AT STEP "<<iterationCounter<<" = \n"<<actSol.ToString()<<std::endl;

      // That's a bit dirty but it's currently the only possible way I see
      algsys_->InitSol(solVecMH_);

      if (IS_LOG_ENABLED(stdsolvestep, dbg3)) {
        ftRes.SetFrequencyResult(actSol);
        ftRes.FourierToTime();
        for(UInt i = 0; i < ftRes.GetNumTimeSteps(); ++i) {
          std::cout<<"TIME RESULT "<<i<<" = "<<ftRes.GetTimeResult(i).ToString()<<std::endl;
        }
      }

      // Calculation relative residual error
      Double residualErr;
      if ( RhsLinL2Norm > 1.0 )
        residualErr = residualL2Norm / RhsLinL2Norm;
      else
        residualErr = residualL2Norm;

      // calculate incremental error
      Double incrementalErr;
      Double solIncrL2Norm = solInc.NormL2();
      Double actSolL2Norm  = actSol.NormL2();

      if ( actSolL2Norm ) incrementalErr = solIncrL2Norm / actSolL2Norm;
      else {
        incrementalErr = solIncrL2Norm;
        WARN("Zero solution vector!! ");
      }

      // Output of norms and data
      if ( nonLinLogging_ == true ) {
        //UInt actStep = PDE_.GetSolveStep()->GetActStep();
        unsigned int solStep = PDE_.IsIterCoupled() ? couplingIter_ : 1;
        OutputNonLinIterInfo(pdename_, solStep, 1, iterationCounter, residualErr, incrementalErr, etaLineSearch);
      }

      // boolean variable, holds condition if another iteration step is necessary
      std::cout << "========= Iterationstep = "<<iterationCounter << std::endl;
      performOneMoreStep = ((incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_)) && (residualErr <= lastError_);
      if(residualErr >= lastError_){
        std::cout<<"Stopped due to divergence at iteration "<< iterationCounter << std::endl;
      }
      lastError_ = residualErr;
      std::cout<<"========= incrementalErr = "<<incrementalErr<<std::endl;
      std::cout<<"========= residualErr = "<<residualErr<<std::endl;
      if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
        EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
            << "' in step no '" << 1
            << "' at iteration '" << iterationCounter
            << "'.\n ==> incremental error: " << incrementalErr
            << "\n ==> residual error: " << residualErr);
      }
    }while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

  }

  void StdSolveStep::EvaluateNonlinearity(MHTimeFreqResult& ftRes,
                                          const SBM_Vector& actSol){

    // Register the multiharmonic solution at MHTimeFreqResult
    ftRes.SetFrequencyResult(actSol);

    // and transform the solution into time-domain to be able
    // to evaluate the nonlinearity (e.g. BH curve in electromagnetics)
    ftRes.FourierToTime();


    // ============================================================================
    // Evaluation of nonlinearity via callback mechanism in CoefFunctionHarmBalance
    // ============================================================================

    // Now we have to evaluate the nonlinearity at every integration point
    // (not yet possible, so we evaluate it element wise)
    // and transform the time signal back into frequency domain.
    // Therefore we loop over every time step, which MHTimeFreqResult provides
    // and set the solVec_ pointer to the correct sub-SBM vector in FeFunctions.
    // Then the CoefFunctionHarmBalance should evaluate the nonlinearity at the
    // integration points and store the solution...still in time domain.
    // The transformation back into the frequency domain is carried out after this
    // loop via changing the math parser variable "finishCash", which is a callback
    // to CoefFunctionHarmBalance::UpdateSolution() and transforms the nu(t) into
    // nu(harmonic)

    // Loop over time steps
    for(UInt i = 0; i < ftRes.GetNumTimeSteps(); ++i){
//      std::cout << "=========== Timestep " << i << " =========" << std::endl;
      // TODO this should be double
      Vector<Complex> timeStepVec = ftRes.GetTimeResult(i);
      solVec_(0) = (Vector<Complex>)timeStepVec;

      // Trigger the callback mechanism in CoefFunctionHarmBalance
      // Now the PDE has the solution vector via the FeSpace and we activate
      // the callback mechanism to cache the solution vector for current harmonic
      mParser_->SetValue(MathParser::GLOB_HANDLER, "cacheResult", i);
    }

    // Now that the nu(t) results are cached, we can perform the FFT
    UInt f = 1;
    mParser_->SetValue(MathParser::GLOB_HANDLER, "finishCash", f);
  }

  void StdSolveStep::AssembleMH(const UInt& N, const UInt& M, const bool onlyDiagBlocks) {
    // loop over every frequency and assemble the correct SBM blocks

//    std::cout << "  - Calculating BiLinearForms for multiharmonic analysis" <<std::endl;

    // Init all matrices, which have to be reassembled
    // Usually this is done in Assemble::AssembleMatrices_Std but we don't
    // use this method, therefore we call a special method here.
    assemble_->InitMultHarm();


    // Special treatment is needed for the diagonal blocks, due to the mass part,
    // therefore handle this case seperately
    // Contrary to the usual assembling process, we have to pass regionNonLinTypes_
    // because regions without a BH curve don't have to be assembled into off-diagonal
    // blocks in the global system matrix...performance improvement
    assemble_->AssembleMatrices_MultHarm(0, solStrat_->GetNumHarmN(),
        solStrat_->GetNumHarmM(),
        regionNonLinTypes_,
        multHarmFreqVec_);


    if(!onlyDiagBlocks){
      for (UInt i = 0; i < multHarmFreqVec_.GetSize(); ++i) {
        Integer tmpH = domain->GetDriver()->HarmonicOfIndex(i);
        Integer h;
        if(domain->GetDriver()->IsFullSystem()){
          // Not optimized version including all harmonics (even and odd ones)
          h = tmpH;

          // the matrix entries for odd harmonics are zero, therefore we don't
          // have to assemble them
          // Have to assemble them in certain cases...
//           if( h%2 != 0 && h!=0 ){
//             continue;
//           }


        }else{
          // Ok, now it gets confusing because in the performance-optimized
          // version, we draw a border between the harmonics of the system matrix
          // and the solution vector
          if(tmpH == 0) h = 0;
          else if(tmpH < 0) h = tmpH - 1;
          else h = tmpH + 1;
        }

        mParser_->SetValue(MathParser::GLOB_HANDLER, "harmonicHandle", h);

        // And set the corresponding frequency
        Double tmpf = (Double)solStrat_->GetBaseFreq();
        Double freq = tmpf * h;
        mParser_->SetValue(MathParser::GLOB_HANDLER, "f", freq);

        if (IS_LOG_ENABLED(stdsolvestep, dbg3)) {
          std::cout<<"harmonic = "<<h<<", frequency = "<<freq<<std::endl;
        }


        if( std::abs(h) > (Integer)M || h == 0) {
          continue;
        } else {
          // Assemble the correct SBM-block, therefore pass the harmonic (-N,...,0,...,N)
          // Contrary to the usual assembling process, we have to pass regionNonLinTypes_
          // because regions without a BH curve don't have to be assembled into off-diagonal
          // blocks in the global system matrix...performance improvement
          // NOTE: In the new optimized version, only odd harmonics are considered
          assemble_->AssembleMatrices_MultHarm(h, solStrat_->GetNumHarmN(), solStrat_->GetNumHarmM(), regionNonLinTypes_);
        }
      }
    }

    // Flag the the matrices were assembled at least once
    mParser_->SetValue(MathParser::GLOB_HANDLER, "harmonicHandle", 0);
  }

  void StdSolveStep::GetSolutionValMultHarm(const UInt& h){
    algsys_->GetSolutionVal(h, solVec_);
  }


  void StdSolveStep::GetRHSValMultHarm(const UInt& h){
    algsys_->GetRHSVal(h, rhsVec_);
  }



  // ======================================================
  // METHODS FOR EIGENVALUE COMPUTATION
  // ======================================================
  
  UInt StdSolveStep::CalcEigenFrequencies( Vector<Double>& frequencies, Vector<Double>& errBounds,
          UInt numFreq, double shift, bool sort) {
    
    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();
    
    assemble_->AssembleMatrices();
    
    // Setup solver
    algsys_->SetupEigenSolver(numFreq, shift, false, sort, false);
    
    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies(frequencies, errBounds);
    
    return frequencies.GetSize();
  }
  
  UInt StdSolveStep::CalcEigenFrequencies( Vector<Complex>& frequencies, Vector<Double>& errBounds,
          UInt numFreq, Double shift, bool sort, bool bloch) {
    
    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();
    
    assemble_->AssembleMatrices();
    
    // Setup solver  - we cannot be quadratic and bloch concurrently!
    algsys_->SetupEigenSolver(numFreq, shift, !bloch, sort, bloch);
    
    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies( frequencies, errBounds );
    algsys_->ExportLinSys(false, false, true); // setup
    
    return frequencies.GetSize();
  }
  
  UInt StdSolveStep::CalcEigenFrequencies( Vector<Double>& frequencies, Vector<Double>& errBounds, Double minVal, Double maxVal){
      algsys_->InitRHS();
      algsys_->InitSol();
      algsys_->InitMatrix();
      assemble_->AssembleMatrices();
      // Setup solver  - we cannot be quadratic and bloch concurrently!
      //algsys_->SetupEigenSolver(minVal, maxVal);
      // Calculate eigenfrequencies
      algsys_->CalcEigenFrequencies( frequencies, errBounds );
      algsys_->ExportLinSys(false, false, true); // setup
      return frequencies.GetSize();
  }

  void StdSolveStep::CalcEigenValues(BaseVector &sol, BaseVector &err, Double minVal, Double maxVal ){
      algsys_->CalcEigenValues( sol, err, minVal, maxVal );
  }

  void StdSolveStep::GetEigenMode( UInt numMode, bool right ) {
    algsys_->GetEigenMode( numMode, right );
    
    // Get the solution and store it
    // Since the entries of solVec_ are pointers to the SingleVector
    // of the FE function, it automatically inserts the values there
    algsys_->GetSolutionVal(solVec_);
  }
  
  
  
  
  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================

  Double StdSolveStep::SetLinRHS( Double loadFactor, bool nonlin, bool multiharmonic)
  {
    
    //std::cout << "SetLinRHS with bool nonlin = " << nonlin << std::endl;
    Double RhsLinL2Norm;
    
    // to incorporate loads
    if(nonlin){
      assemble_->AssembleNonLinRHS();
    } else {
      assemble_->AssembleLinRHS();
    }
    
    //Set special RHS Values
    PDE_.SetRhsValues();
    
    // Stores rhs vector into extForces and returns that L2-norm
    if(multiharmonic) algsys_->GetFullMultiHarmRHSVal(RhsLinVal_);
    else algsys_->GetRHSVal( RhsLinVal_ );
    
    RhsLinVal_.ScalarMult(loadFactor);
    
    RhsLinL2Norm = RhsLinVal_.NormL2();
    
    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      // Note: there are PDEs, such as elecconduction, which always have rhs=0. Those should not emit this warning. SE.
      if (pdename_ != "elecConduction" && pdename_ != "electrostatic"){
        WARN("Zero external force vector!! ");
      }
    }
    
    return RhsLinL2Norm;
  }
  
  UInt StdSolveStep::SetDeltaLinRHS()
  {
    
    // to incorporate loads
    assemble_->AssembleLinRHS(); 
    //Set special RHS Values
    PDE_.SetRhsValues();
    
    SBM_Vector newRhsLinVal(BaseMatrix::DOUBLE); //!< external forces (for nonlin simulations)
    algsys_->GetRHSVal( newRhsLinVal );
    DeltaRhsLinVal_.Add( 1.0, newRhsLinVal,
            -1.0, tmpOldRhsLinVal_ );
    
    RhsLinVal_=tmpOldRhsLinVal_;
    
    Double DeltaNorm=DeltaRhsLinVal_.NormL2();
    Double oldNorm  =tmpOldRhsLinVal_.NormL2();
    Double aux;
    UInt nrLoadSteps;
    
    UInt minNrLoadSteps = 1;
    UInt maxNrLoadSteps = 1;
    
    if(oldNorm<1e-13)
      nrLoadSteps=maxNrLoadSteps;
    else{
      aux= Double(maxNrLoadSteps)*(DeltaNorm/oldNorm);
      nrLoadSteps=UInt(aux);
    }
    if(nrLoadSteps<minNrLoadSteps)
      nrLoadSteps=minNrLoadSteps;
    else if (nrLoadSteps>maxNrLoadSteps)
      nrLoadSteps=maxNrLoadSteps;
    
    oldRhsLinVal_=tmpOldRhsLinVal_;
    tmpOldRhsLinVal_=newRhsLinVal;
    
    return nrLoadSteps;
  }
  
  Double StdSolveStep::LineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol,
          Double& etaLineSearch, bool trans)  {
    
    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 4;
    const Double eta[nrEtas] = {0.1, 0.25, 0.5, 1.0}; //, 0.5, 0.25, 0.125, 0.1}; // TODO check why we search from small to large?!
    
    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;
    
    for( UInt i=0; i<nrEtas; i++) {
      //std::cout << "Testing eta = " << eta[i] << std::endl;
      // take care about the data types
      if(actSol.GetEntryType() == BaseMatrix::DOUBLE) actSol.Add( 1.0, solOld, eta[i], solIncrement);
      else actSol.Add( (Complex) 1.0, solOld, (Complex) eta[i], solIncrement);

      //store new solution
      solVec_ = actSol;

      // set RHS: linear part
      algsys_->InitRHS(RhsLinVal_ );
      // and nonlinpart if any
      assemble_->AssembleNonLinRHS();


      // setup the matrices
      bool isNewton = false;
      assemble_->AssembleMatrices(isNewton);


      if( trans ) {
        //now update RHS according to time stepping
        std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
        std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
        std::map<FEMatrixType,Integer>::iterator matIt;
        UInt pos = 0;
        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
          if(matIt->second < 0)
            continue;
          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
            fncIt->second->GetTimeScheme()->ComputeStageRHS(0,matIt->second,stageRHS_.GetPointer(pos));
          }
          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
        }

        //substract from RHS the term K*sol
        solVec_.ScalarMult(-1.0);
        algsys_->UpdateRHS(STIFFNESS,solVec_,true);
        algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
        solVec_.ScalarMult(-1.0);
      }
      else {
        solVec_.ScalarMult(-1.0);
        algsys_->UpdateRHS(SYSTEM,solVec_,true);
        solVec_.ScalarMult(-1.0);
      }


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      SBM_Vector actRHS(BaseMatrix::DOUBLE);
      algsys_->GetRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = actRHS.NormL2();

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
      LOG_DBG(stdsolvestep) << "LS eta=" << eta[i] << " norm=" << residualL2Norm << " eta_opt=" << etaOpt << " norm_opt=" << residualL2NormOpt;
    }

    //std::cout << "Optimal eta = " << etaOpt << std::endl;
    etaLineSearch = etaOpt;

    // Set new solution
    actSol.Add( 1.0, solOld, etaOpt, solIncrement );

    return residualL2NormOpt;
  }


  Double StdSolveStep::LineSearchMultHarm(const SBM_Vector& solIncrement, SBM_Vector& actSol,
          Double& etaLineSearch, MHTimeFreqResult& ftRes)  {

    SBM_Vector solOld(BaseMatrix::COMPLEX);
    solOld = actSol;

    UInt h;
    h = solStrat_->GetNumHarmN();

    const UInt nrEtas = 4;
    const Double eta[nrEtas] = {0.1, 0.25, 0.5, 1}; //, 0.5, 0.25, 0.125, 0.1};

    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    for( UInt i=0; i<nrEtas; i++) {

      LOG_DBG(stdsolvestep) <<" LineSearchMultHarm: Testing eta = " << eta[i];


      if(actSol.GetEntryType() == BaseMatrix::DOUBLE){
        EXCEPTION("StdSolveStep::LineSearchMultHarm Solution vector is real valued in multiharmonic analysis!");
      }

      if(lineSearch_ == "multiharmonicIncreasing"){
        //this is the adv linesearch - this doesnt improve the speed but the quality of the result
        // add only the considered harmonics while linesearching
        for(UInt x=0; x<=consideredH_;x++){
          if(x==0){
            actSol.Add( (Complex) 1.0, solOld, (Complex) eta[i], solIncrement,h);
          }else{
            actSol.Add( (Complex) 1.0, solOld, (Complex) eta[i], solIncrement,h+x);
            actSol.Add( (Complex) 1.0, solOld, (Complex) eta[i], solIncrement,h-x);
          }
        }
      } else { // standard linesearch
          // Actually it's this actSol = solOld + eta[i] * solIncrement;
          actSol.Add( (Complex) 1.0, solOld, (Complex) eta[i], solIncrement);
      }

      // We need to do this in order to evaluate at the correct result vector
      algsys_->InitSol(actSol);

      // Evaluate the nonlinearity (e.g. BH curve in electromagnetics)
      this->EvaluateNonlinearity(ftRes, actSol);

      // set RHS: linear part
      algsys_->InitRHS(RhsLinVal_ );
      // and nonlinpart if any
      //assemble_->AssembleLinRHS();


      // setup the matrices
      this->AssembleMH(solStrat_->GetNumHarmN(), solStrat_->GetNumHarmM());
      assemble_->PostAssemble();

      // Computation of effective matrix:
      /* NOTE: this is commented because we also include the MASS matrix
               in the SYSTEM matrix, as defined in Assemble::CreateMatrixMap().
               Sometimes having an extra MASS matrix is benefitial, e.g. for exporting and
               comparing different matrix parts, that's why it's still here
      */
      //std::map<FEMatrixType,Double> empty;
      //algsys_->ConstructEffectiveMatrix(NO_FCT_ID,  empty, true );



      actSol.ScalarMult(-1.0);
      algsys_->UpdateRHS_MultHarm(SYSTEM,actSol,true);
      actSol.ScalarMult(-1.0);

      // =====================================================================
      // calculation of error norms
      // =====================================================================
      SBM_Vector actRHS(BaseMatrix::COMPLEX);
      algsys_->GetFullMultiHarmRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = actRHS.NormL2();

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }else if(residualL2Norm > residualL2NormOpt){
        continue;
        //TODO Proof Theory that  the first lokal minimum is the global minimum of all residuals as well.
      }
    } // eta-loop

//    std::cout << "Optimal eta = " << etaOpt << std::endl;
    etaLineSearch = etaOpt;

//    std::cout << "Considered Harmonics: " << consideredH_ << std::endl;
//    std::cout << eta[0] << " - " << etaError[0] << std::endl;
//    std::cout << eta[1] << " - " << etaError[1] << std::endl;
//    std::cout << eta[2] << " - " << etaError[2] << std::endl;
//    std::cout << eta[3] << " - " << etaError[3] << std::endl;

    // Set new solution for optimal eta
    if(lineSearch_ == "multiharmonicIncreasing"){
      if(residualL2NormOpt > lastError_){
        //if error is bigger than error from last time,
        // increase the considered harmonics
        // And do not change the residuum -> set etaOpt to zero
        if(consideredH_ < h){
          consideredH_++;
        }
        etaOpt = 0;
        //update error
        residualL2NormOpt = lastError_;
        //we do not have to update the solution
        return residualL2NormOpt;
      } else {
        etaLineSearch = etaOpt;
        for(UInt x=0; x<=consideredH_ ;x++){
          if(x==0){
            actSol.Add( (Complex) 1.0, solOld, (Complex) etaOpt, solIncrement,h);
          }else{
            actSol.Add( (Complex) 1.0, solOld, (Complex) etaOpt, solIncrement,h+x);
            actSol.Add( (Complex) 1.0, solOld, (Complex) etaOpt, solIncrement,h-x);
          }
        }
      }
    } else{ //standard linesearch
      etaLineSearch = etaOpt;
      actSol.Add( (Complex)1.0, solOld, etaOpt, solIncrement );
    }
    return residualL2NormOpt;
  }
  
  

  Double StdSolveStep::LineSearchMag(SBM_Vector& solIncrement, SBM_Vector& actSol,
          Double& etaLineSearch, bool trans)
  {
    
    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 4;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.1};
    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;
    
    for( UInt i=0; i<nrEtas; i++) {
      actSol.Add( 1.0, solOld, eta[i], solIncrement);
      //      actSol = solIncrement * eta[i];
      //      actSol += solOld;
      
      //store new solution
      solVec_ = actSol;
      
      // recalculate RHS with new values to get new residual (f^(k+1))========
      algsys_->InitRHS(RhsLinVal_ );
      
      if( trans ) {
        //EXCEPTION("Line Search for nonlinear transient problems yet verified")
        //        assemble_->AssembleNonLinRHS();
        //        TS_alg_->UpdateRHS(actSol);
        assemble_->AssembleNonLinRHS();
        //PDE_.SetRhsValues();
        //now update RHS according to time stepping
        std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
        std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
        std::map<FEMatrixType,Integer>::iterator matIt;
        UInt pos = 0;
        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
          if(matIt->second < 0)
            continue;
          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
            fncIt->second->GetTimeScheme()->ComputeStageRHS(0,matIt->second,stageRHS_.GetPointer(pos));
          }
          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
        }
      }
      else {
        assemble_->AssembleNonLinRHS();
      }
      
      
      // =====================================================================
      // calculation of error norms
      // =====================================================================
      SBM_Vector actRHS(BaseMatrix::DOUBLE);
      algsys_->GetRHSVal( actRHS );
      
      // calculation of residual error =======================================
      Double residualL2Norm = actRHS.NormL2(); // L2Norm of  ( f_i^(k+1) - f_a )
      
      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }
    
    etaLineSearch = etaOpt;
    
    // Careful: in the end, we have to re-assemble the RHS with the correct
    // value i.e. use the "optimal" solution
    actSol.Add(1.0, solOld, etaOpt, solIncrement );
    
    return residualL2NormOpt;
  }
  
  Double StdSolveStep::LineSearchMaterial(SBM_Vector& solIncrement, SBM_Vector& actSol,
          Double& etaLineSearch, Double& RhsLinL2Norm, bool trans)
  {
    REFACTOR;
    /*
     SBM_Vector solOld;
     solOld = actSol;
     const UInt nrEtas = 3;
     const Double eta[nrEtas] = {0.9, 0.5, 0.3};
     //    const Double eta[nrEtas] = {0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1.0};
     // initialize etaOpt or receive compiler warning
     Double etaOpt = 0.0;
     Double residualL2NormOpt = 1e15;
     
     SBM_Vector tmpSol(BaseMatrix::DOUBLE);
     algsys_->GetSolutionVal( tmpSol );
     
     for( UInt i=0; i<nrEtas; i++) {
     //       if (i>0)
     //         actSol=solOld;
     actSol.Add( 1.0, solOld,
     eta[i], solIncrement );
     //actSol = solOld + solIncrement * eta[i];
     //      actSol -= solOld;
     
     
     //store new solution
     tmpSol = actSol;
     solVec_ = tmpSol;
     
     // Recalculate residual, f-Cu-Mu-K*u
     algsys_->InitRHS();
     
     assemble_->AssembleLinRHS();
     
     // assemble!
     assemble_->AssembleMatrices();
     
     // account for Dirichlet BCs
     PDE_.SetBCs();
     
     algsys_->ConstructEffectiveMatrix(matrix_factor_);
     
     algsys_->BuildInDirichlet();
     
     TS_alg_->UpdateRHS(actSol);
     // substract K^* u^k from RHS
     
     algsys_->RemoveIDBCInfoFromMatrix();
     
     TS_alg_->SubstractStiffnessFromRHS(actSol);
     
     
     // =====================================================================
     // calculation of error norms
     // =====================================================================
     SBM_Vector actRHS(BaseMatrix::DOUBLE);
     algsys_->GetRHSVal( actRHS );
     
     // calculation of residual error =======================================
     Double residualL2Norm = actRHS.NormL2(); // L2Norm of  (f-Ku )
     
     Double residualErr;
     
     if ( RhsLinL2Norm > 1.0 )
     residualErr    = residualL2Norm /  RhsLinL2Norm;
     else
     residualErr    = residualL2Norm;
     
     if (residualL2Norm < residualL2NormOpt) {
     residualL2NormOpt = residualL2Norm;
     etaOpt = eta[i];
     }
     }
     
     etaLineSearch = etaOpt;
     
     actSol.Add( 1.0, solOld, etaOpt, solIncrement );
     
     return residualL2NormOpt;
     */
    return 0.0;
  }
  
  // read nonlinear parameters from xml file
  void StdSolveStep::ReadNonLinData()
  {
    // Get ParamNode of pde
    PtrParamNode nonLinNode = solStrat_->GetNonLinNode();
    
    nonLinInfo_ = PDE_.GetInfoNode()->Get("nonlinearConvergence");

    static_non_lin_step_timer_ = nonLinInfo_->Get("timer")->AsTimer();
    static_non_lin_step_timer_->SetSub();

    // Check, if any nonlinear node was found
    if( !nonLinNode ) {
      WARN("Taking default parameters for nonlinear data" );
    }
    
    // Read data, if "nonLinear" element was found
    if(nonLinNode)
    {
      // solution method
      nonLinNode->GetValue( "method", nonLinMethod_, ParamNode::PASS );
      
      // perform logging?
      nonLinNode->GetValue( "logging", nonLinLogging_, ParamNode::PASS );
      if(nonLinLogging_)
      {
        if(logFile_.is_open())
          logFile_.close();
        logFile_.open("nonlin.txt");
        logFile_ << "#1:iterationCounter \t2:residualErr \t3:incrementalErr \t4:etaLineSearch" << std::endl;
      }
      
      // type of line search
      if( nonLinNode->Has("lineSearch") ) {
        nonLinNode->Get( "lineSearch")->GetValue( "type", lineSearch_,ParamNode::PASS );
      }
      
      // incremental stopping criterion
      nonLinNode->GetValue( "incStopCrit", incStopCrit_, ParamNode::PASS );
      
      // residual stopping criterion
      nonLinNode->GetValue( "resStopCrit", residualStopCrit_, ParamNode::PASS );
      
      // maximal number of NL-iterations
      nonLinNode->GetValue( "maxNumIters", nonLinMaxIter_, ParamNode::PASS );
      
      // abort if max number of iterations is reached?
      nonLinNode->GetValue("abortOnMaxNumIters",abortOnMaxIter_,ParamNode::INSERT);
      
      LOG_TRACE(stdsolvestep) << "Nonlinear convergence criteria were read:";
      LOG_DBG3(stdsolvestep) << "\tincremental Stopping Criterion: " << incStopCrit_;
      LOG_DBG3(stdsolvestep) << "\tresidual Stopping Criterion: " << residualStopCrit_;
      LOG_DBG3(stdsolvestep) << "\tmaxNumIters: " << nonLinMaxIter_;
      if( nonLinNode->Has("stopOnLimit") ) {
        std::string solutionString;
        // quantity
        solutionString = nonLinNode->Get( "stopOnLimit")->Get( "quantity" )->As<std::string>();
        solutionLimit_ = SolutionTypeEnum.Parse(solutionString);
        if (feFunctions_.find(solutionLimit_) == feFunctions_.end() ) 
          EXCEPTION("ERROR: Solution type '" << solutionString << "' is not part of PDE '" << pdename_ << 
                  "' and cannot serve as stopping limit criterion");
        // minimum value
        nonLinNode->Get( "stopOnLimit")->GetValue( "min", minValidValue_, ParamNode::PASS );
        // maximum value
        nonLinNode->Get( "stopOnLimit")->GetValue( "max", maxValidValue_, ParamNode::PASS );
        // region
        std::string solutionRegion;
        solutionRegion = nonLinNode->Get( "stopOnLimit")->Get( "region" )->As<std::string>();
        solutionLimitReg_ = ptgrid_->GetRegion().Parse(solutionRegion);
        if (solutionLimitReg_ != ALL_REGIONS) {
          if( subdoms_.Find(solutionLimitReg_) == -1 ) 
            EXCEPTION("ERROR: Region '" << solutionRegion <<"' is not part of PDE '" << pdename_ << 
                    "' and cannot serve as stopping limit region");
          
          // don't know how to implement region checking as we only have the solution as vector (glmvector_[0])
          EXCEPTION("checking limits in a region is not implemented (and I don't know how) (SE). " 
                  << "Just remove the 'region=XXX' tag");
        }
        LOG_DBG3(stdsolvestep) << "\tStop on Limit: " << solutionString;
        LOG_DBG3(stdsolvestep) << "\t\tmin:     " << minValidValue_;
        LOG_DBG3(stdsolvestep) << "\t\tmax:     " << maxValidValue_;
        LOG_DBG3(stdsolvestep) << "\t\tin Region:     " << solutionRegion;
        LOG_DBG3(stdsolvestep) << "\t\tRegionID:     " << solutionLimitReg_;
      }
    }

    PtrParamNode setup = nonLinInfo_->Get("setup");
    setup->Get("method")->SetValue(nonLinMethod_);
    setup->Get("lineSearch")->SetValue(lineSearch_);
    setup->Get("maxIter")->SetValue(nonLinMaxIter_);
    setup->Get("incStopCrit")->SetValue(incStopCrit_);
    setup->Get("residualStopCrit")->SetValue(residualStopCrit_);
    setup->Get("logging")->SetValue(nonLinLogging_);
  }
  
  void StdSolveStep::OutputNonLinIterInfo(const std::string& pdeName,
          UInt solStep, UInt iterationCounter, Double residualErr,
          Double incrementalErr, double etaLineSearch, int coupledIterStep)
  {
    assert(nonLinInfo_);
    PtrParamNode nlc = nonLinInfo_;

    // usually we have only one solStep.
    // In the detail case we have repeatet solStep value="1" with different analysis_id
    // In the non-detail case we overwrite the analysis_id but also need to clear the old content.
    string aid = domain->GetDriver()->GetAnalysisId().ToString();
    PtrParamNode ss;

    if(progOpts->DoDetailedInfo())
      ss = nlc->GetByVal("solStep", "value", std::to_string(solStep), "analysis", aid);
    else
    {
      ss = nlc->GetByVal("solStep", "value", std::to_string(solStep));
      ss->Get("analysis")->SetValue(aid); // possibly overwrite
    }

    // in the coupling step we have a layer between solStep and iteration
    if(coupledIterStep >= 0)
      ss = ss->GetByVal("couplingStep", "value", coupledIterStep, ParamNode::INSERT);

    if(!progOpts->DoDetailedInfo() && iterationCounter <= 1)
      ss->ClearChildren("iteration"); // do it here to delete only iterations and not other coupling steps

    PtrParamNode iter = ss->Get("iteration",ParamNode::APPEND);
    iter->Get("pdeName")->SetValue(pdeName);
    iter->Get("nr")->SetValue(iterationCounter);
    if(residualErr >= 0)
      iter->Get("residualErr")->SetValue(residualErr);
    if(incrementalErr >= 0)
      iter->Get("incrementalErr")->SetValue(incrementalErr);
    if(etaLineSearch > 0)
      iter->Get("eta_linesearch")->SetValue(etaLineSearch);

    if(nonLinLogging_)
    {
      logFile_ << iterationCounter << "\t"
               << residualErr << "\t"
               << incrementalErr << "\t"
               << etaLineSearch << std::endl;
      logFile_.flush();
    }

  }

  bool StdSolveStep::GetRidOfZerosActive(){
    // Note: Currently we do not support multi-grid or static condensation with this approach

    // default parameters  (false by default)
    bool supportedBySolver = true;
    string getRidOfZerosXML;
    bool hasNCI = false;
    
    PtrParamNode myParam = this->PDE_.GetDomain()->GetParamRoot();
    PtrParamNode seqStepParamNode = myParam->GetByVal("sequenceStep", std::string("index"),
                          this->PDE_.GetDomain()->GetDriver()->GetActSequenceStep());

    // get info if we have NCI or not
    if (this->ptgrid_->HasNCI()) {
        hasNCI = true;
    }
    
    // get XML input ("auto" and 1e-20 are defaluts from XML)
    getRidOfZerosXML = seqStepParamNode->Get("linearSystems/getRidOfZeros")->As<string>();
    getRidOfZerosTol_ = seqStepParamNode->Get("linearSystems")->Get("getRidOfZerosTolerance")->As<Double>();

    // now check the solver list if we use a solver that is supported
    if ( algsys_->GetSolver()->GetSolverType() != BaseSolver::PARDISO_SOLVER ) {
      supportedBySolver = false;
    }
    
    ParamNodeList analysisNode = seqStepParamNode->Get("analysis")->GetChildren();
    std::string analysisName = analysisNode[0]->GetName();

    // now adapt the switch if any of the untested / not working cases
    std::string useCase = "";
    if( algsys_->UseAMG() )
      useCase = "AMG";
    else if ( algsys_->UseStaticCondensation() )
      useCase = "static condensation";
    else if ( !supportedBySolver || algsys_->HasPrecond() )
      useCase = "solver setup";
    else if ( analysisName == "inverseSource" )
      useCase = "inverseSource";
    else if ( myParam->Has("optimization") )
      useCase = "optimization";
    else if (this->PDE_.IsNonLin())
      useCase = "nonLinear";

    // we select the scenario based on gatherd data
    if (useCase != "" && hasNCI && getRidOfZerosXML=="auto") {
      std::cout << "StdSOlveStep::GetRidOfZeros: This feature is not available for the current use case (" << useCase << ")"; 
      return false;
    }
    else if (useCase == "" && hasNCI && getRidOfZerosXML=="auto") {
      std::cout << "Zero entities will be removed from the system matrix in each iteration to reduce solver effort because the model contains at least one NCI. Define \"no\" for \"getRidOfZeros\" in \"linearSystems\" explicitly to avoid this.";
      return true;
    }
    else if (useCase != ""  && getRidOfZerosXML=="yes") {
      Exception("StdSOlveStep::GetRidOfZeros: This feature is not available for the current use case (" + useCase + ")");
      return false;
    }
    else if (useCase == "" && getRidOfZerosXML=="yes") {
      return true;
    }
    else if (getRidOfZerosXML=="no") {
      return false;
    }
    else {
      return false;
    }
    
  }
  
} // end of namespace
