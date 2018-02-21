#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>

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
#include "DataInOut/Logging/LogConfigurator.hh"
//#include "MatVec/SingleVector.hh"


namespace CoupledField {
  
  // declare logging stream
  DECLARE_LOG(stdsolvestep)
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
    
    // Copy vectors FE functions in SBM-vector for communication
    // with OLAS and time stepping
    solVec_.SetSize( feFunctions_.size() );
    rhsVec_.SetSize( feFunctions_.size() );
    
    //    resVec_.SetSize( feFunctions_.size() );
    //    nonLinRHS_.SetSize( feFunctions_.size() );
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator it;
    it = feFunctions_.begin();
    // UInt pos = 0;
    for( ; it != feFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      FeFctIdType id = ptFct->GetFctId();
      solVec_.SetSubVector(ptFct->GetSingleVector(), id);
    }
    //pos = 0;
    it = rhsFeFunctions_.begin();
    for( ; it != rhsFeFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      FeFctIdType id = ptFct->GetFctId();
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
    
    //		std::cout << "StdSolveStep - Constructor: " << std::endl;
    //		std::cout << "pdename_: " << pdename_ << std::endl;
    //		std::cout << "isHyst_: " << isHyst_ << std::endl;
    
    /*
     * for hysteresis, initialize additional values needed by
     * CalcResidualAndConfigSystemForHysteresis
     */
    if(isHyst_){
      
      currentLinRhsVec_ = rhsVec_;
      // set all rhs and res values to rhsVec_ first to initialize the
      // correct size
      currentNonLinRhsVec_ = rhsVec_;
      currentNonLinRhsVec_.Init();
      
      currentRHSload_partial_ = rhsVec_;
      currentRHSload_partial_.Init();
      
      currentResVec_ = currentNonLinRhsVec_;
      
      oldTSLinRhsVec_ = currentNonLinRhsVec_;
      oldTSNonLinRhsVec_ = currentNonLinRhsVec_;
      
      oldItResVec_ = currentNonLinRhsVec_;
      oldItNonLinRhsVec_ = currentNonLinRhsVec_;
      
      oldTSSolVec_ = solVec_;
      oldTSSolVec_.Init();
    }
    
    // In the end, read nonlinear data from xml-file
    if( nonLin_ || nonLinMaterial_ ) {
      ReadNonLinData();
    }
    
    //logFile_.open("nonlin.txt");
    
    mHandle_ = PDE_.GetDomain()->GetMathParser()->GetNewHandle();
    mParser_ = PDE_.GetDomain()->GetMathParser();
    mParser_->SetExpr(mHandle_,"step");
  }
  
  
  //! Destructor
  StdSolveStep::~StdSolveStep() {
    //logFile_.close();
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
      algsys_->SetupPrecond();
      
      algsys_->SetupSolver();
    }
    
    // Solve problem
    algsys_->Solve();
    
    // Get the solution and store it
    algsys_->GetSolutionVal(solVec_);
  }
  
  
  void StdSolveStep::StepStaticNonLin() {
    
    bool performOneMoreStep;
    bool isNewton = false;
    
    SBM_Vector solInc(BaseMatrix::DOUBLE);
    
    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;
    
    // =================================
    //  Outer loop: Multilevel strategy 
    // =================================
    UInt numLevels = solStrat_->GetNumSolSteps();
    for( UInt iLevel = 0; iLevel < numLevels; ++iLevel ) {
      
      // create new timer object and put it to related info element
      shared_ptr<Timer> timer(new Timer());
      PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
      iter->GetByVal("solStep","value",iLevel+1,ParamNode::INSERT)
      ->Get("timer")->SetValue(timer);
      timer->Start();
      
      // update the current solution step in a multilevel approach and
      // inform PDEs (containing the FeSpaces), as well as the AlgebraicSystem
      solStrat_->SetActSolStep(iLevel + 1);
      ReadNonLinData();
      PDE_.UpdateToSolStrategy();
      algsys_->UpdateToSolStrategy();
      
      // set the boundary conditions
      PDE_.SetBCs();
      
      //perform the load-steps
      Double loadFactor = 1.0;
      PDE_.GetInfoNode()->Get("PDE")->Get(pdename_)->
      Get("load_factor")->SetValue(loadFactor);
      
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
        
        if ( lineSearch_ != "none" || iterationCounter == 1) {
          //add linear right hand side
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
        }
        
        
        // assemble Newton bilinear forms
        isNewton = true;
        assemble_->AssembleMatrices(isNewton);
        
        //compute effective matrix
        algsys_->ConstructEffectiveMatrix( NO_FCT_ID,
                matrix_factor_[NO_FCT_ID] );
        
        algsys_->BuildInDirichlet();
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        
        bool setIDBC = false;
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          setIDBC = true;
        
        algsys_->Solve(setIDBC);
        
        // new solution is only an increment of the full solution =============
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
        
        // output of norms and data
        if ( nonLinLogging_ == true ) {
          //UInt actStep = PDE_.GetSolveStep()->GetActStep();
          
          if (PDE_.IsIterCoupled()) {
            WriteNonLinIterToInfoXML(pdename_, couplingIter_, iLevel+1, iterationCounter, residualErr, incrementalErr, etaLineSearch);
          } else {
            WriteNonLinIterToInfoXML(pdename_, iLevel+1, iterationCounter, residualErr, incrementalErr, etaLineSearch);
          }
          
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
                  << residualErr << "\t"
                  << incrementalErr << "\t"
                  << etaLineSearch << std::endl;
        }
        
        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
                (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
        
        
        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                  << "' in step no '" << iLevel+1 
                  << "' at iteration '" << iterationCounter 
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
      
      // stop timer
      timer->Stop();
    } // loop over levels
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
      StepTransNonLinHysteresis();
    }
    
    //currently not supported
    //	else if ( nonLin_ && nonLinMaterial_ ) {
    //      StepTransNonLinMaterial();
    //    }
    // do a nonlinear time step
    else if (nonLin_ ){//|| PDE_.IsHysteresis_Fixpoint() == true){
      if ( nonLinTotalFormulation_ ){ //|| PDE_.IsHysteresis_Fixpoint() == true){
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
  
  
  void StdSolveStep::StepTransLin()
  {
    //TODO: add consistency check here
    //basically loop over all functions and check if the solution order is the same...
    
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    
    std::map<FEMatrixType,Integer>::iterator matIt;
    
    bool effectiveMatrixUpdated = false;
    
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 ); 
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }
    
    for(UInt i=0;i<numStages;i++){
      effectiveMatrixUpdated = false;
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
      
      
      // if we want to use static condensation we have to perform the timestepping on element level
      if(algsys_->UseStaticCondensation()){
        
        matrix_factor_.clear();
        // get timeintegration factors and store them in a map
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()
          ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
        }
        // assemble all matrix parts (needed to update rhs) and also calculate complete
        // system matrix
        assemble_->AssembleMatrices_CondTrans(false,i,matrix_factor_);
        
        if(assemble_->IsMatrixUpdated()){
          
          //std::cout << "in SolveStepLin: new matrices computed" << std::endl;
          
          // the system matrix was already created so do not init it!
          //        algsys_->InitMatrix(SYSTEM);
          
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
          //        matrix_factor_[fctId]
          effectiveMatrixUpdated = true;
        }
        
      } else {
        
        
        assemble_->AssembleMatrices();
        if(assemble_->IsMatrixUpdated()){
          //if AMG is used, rebuild auxiliary matrix
          auxSet_ = false;
          //std::cout << "in SolveStepLin: new matrices computed" << std::endl;
          
          // set system matrix to zero initially, as ConstructEffectiveMatrix only
          // sums up the contributions
          algsys_->InitMatrix(SYSTEM);
          matrix_factor_.clear();
          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
            FeFctIdType fctId = fncIt->second->GetFctId();
            fncIt->second->GetTimeScheme()
            ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
            algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
          }
          effectiveMatrixUpdated = true;
        }
        
      }
      
      //now compute the effective right hand side
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
        if(matIt->second < 0)
          continue;
        
        
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(fctId));
        }
        algsys_->UpdateRHS(matIt->first,stageRHS_,effectiveMatrixUpdated);
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
      
      if( effectiveMatrixUpdated ){
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
      }
      
      algsys_->Solve();
      
      algsys_->GetSolutionVal(stageSol);
    }
    
    //update stage
    
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end(); ++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep(  );
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
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }
    
    for(UInt i=0;i<numStages;i++){
      //do initialization 
      rhsVec_.Init();
      LOG_DBG(stdsolvestep) << "StepTransNonLin: Stage: " << i ;
      
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
      Double RhsLinL2Norm = SetLinRHS(loadFactor);
      
      // set iteration counter
      UInt iterationCounter=0;
      
      do {
        iterationCounter++;
        
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
              fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
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
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        
        // just set inh. Dirichlet BCs for the first iteration
        bool setIDBC = false;
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          setIDBC = true;
        
        algsys_->Solve(setIDBC);
        // if setIDBC is true, solInc will contain the inhom. Dirichlet values
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
              fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
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
        
        // output of norms and data
        if ( nonLinLogging_ == true ) {
          // get current step 
          UInt actStep = PDE_.GetSolveStep()->GetActStep();
          
          if (PDE_.IsIterCoupled()) {
            WriteNonLinIterToInfoXML(pdename_, couplingIter_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
          } else {
            WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
          }
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
                  << residualErr << "\t"
                  << incrementalErr << "\t"
                  << etaLineSearch << std::endl;
        }
        
        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
                (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
        
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
        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
          //const SingleVector * subv = solVec_.GetPointer(pos);
          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(pos)));
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
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
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
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        
        //always set inhomog. Dirichlet BCs
        bool setIDBC = true;
        
        algsys_->Solve(setIDBC);
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
        
        // output of norms and data
        if ( nonLinLogging_ == true ) {
          WriteNonLinIterToInfoXML(pdename_, 1,iterationCounter, residualErr, incrementalErr, etaLineSearch);
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
                  << residualErr << "\t"
                  << incrementalErr << "\t"
                  << etaLineSearch << std::endl;
        }
        
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
  
  void StdSolveStep::SetLastItOrLastTSSBMVectors(bool lastTS){
    /*
     *  Function needed in NonLinHysteresis case
     *
     *  it stores the current values of rhs, residual, solution, ...
     *  for the next iteration or next ts by copying it to the corresponding
     *  SBMVectors
     */
    
    if(lastTS){
      // oldTS > values after last iteration of previous TS
      // to be stored after iteration suceeded
      oldTSLinRhsVec_ = currentLinRhsVec_;
      oldTSNonLinRhsVec_ = currentNonLinRhsVec_;
      oldTSSolVec_ = solVec_;
    } else {
      // oldIt > values of the current TS but from previous It
      // during first iteration of a new TS, these vectors contain the values
      // after the last iteration of the previous TS (similar as oldTS...)
      oldItNonLinRhsVec_ = currentNonLinRhsVec_;
      oldItResVec_ = currentResVec_;
    }
  }
  
//  
//  void StdSolveStep::CalcResidualAndConfigSystemForHysteresis(SBM_Vector& oldSolution,SBM_Vector& solIncrement,SBM_Vector& stageSol, Double usedEta,
//          UInt stage, UInt callingCnt, UInt evalVersion, bool trans, bool forceReevaluation,
//          bool skipReassembly, bool debugOutput, bool reset){
//    /*
//     * new parateter added 21.7: skipReassembly
//     * > if true, residual will be calculated but new system will not get assembled
//     * > useful during line search as we do not use the newly assembled system in that case (we simply overwrite it during
//     * next eta)
//     * > set this flag to false after the final eta has been found
//     *
//     *
//     * reset will have the same effect as callingCnt = 0
//     */
//    
//    //UInt callingCntBackup = callingCnt;
//    
//    
//    if(debugOutput){
//      std::cout << "<<<< CalcResidualAndConfig >>>>" << std::endl;
//      std::cout << "CallingCnt: " << callingCnt << std::endl;
//      std::cout << "evalVersion: " << evalVersion << std::endl;
//    }
//    
//    // IMPORTANT: reallow reevaluation of hyst operator (as input may have changed)
//    PDE_.SetFlagInCoefFncHyst("resetReeval",true);
//    
//    
//    if(reset){
//      callingCnt = 0;
//    }
//    
//    /*
//     * flags to be set in hysteresis operator to enable different
//     * evalVersions
//     */
//    /*
//     * tensorToReturn---
//     * 1: assemble matrix with epsInitial / nuInitial
//     * 2: assemble matrix with eps0 / nu0
//     * 3: assemble deltaMat
//     */
//    UInt tensorToReturnLHS;
//    UInt tensorToReturnRHS;
//    
//    /*
//     * tensorToAdd---
//     * (only needed if tensorToReturn = 3)
//     * 1: add epsInitial to deltaP/deltaE / add nuInitial to deltaM/deltaH
//     * 2: add eps0 to deltaP/deltaE / add nu0 to deltaM/deltaH
//     */
//    UInt tensorToAddLHS;
//    UInt tensorToAddRHS;
//    
//    /*
//     * vectorToReturn
//     * 1: put P(E_k) / M(H_k) on rhs
//     * 2: put P=0 / M = 0 on rhs
//     * 3: put P(E_k) - P(E_lastTS) on rhs / M(H_k) - M(H_lastTS)
//     * 4: put only P(E_lastTS) / M(H_lastTS) on rhs
//     */
//    UInt vectorToReturn;
//    
//    /*
//     * useLastTS
//     * (only needed if tensorToReturn = 3)
//     * true: deltaMat_k = (P(E_k) - P(E_lastTS))/(E_k - E_lastTS)
//     *        or for magnetics with M and H
//     *
//     * false: deltaMat_k = (P(E_k) - P(E_k-1))/(E_k - E_k-1)
//     *         or for magnetics with M and H
//     */
//    bool useLastTS;
//    
//    /*
//     * evaluationPurpose
//     * 1: assemble
//     * 2,3: store
//     * 4: output
//     *
//     * > here of course only option 1 is used
//     * (for more details see coefFncHyst)
//     */
//    UInt evaluationPurpose = 1;
//    PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
//    
//    /*
//     * flag indicating if we need to reassemble system
//     * (only needed if matrixes on rhs and lhs are different)
//     */
//    bool reassemble;
//    
//    /*
//     * flag for easier switching between different rhs versions
//     * 1: full (re-)evaluation
//     *      electrostatics:
//     *        res_k+1 = f_lin + f_nonlin( E_k+1 ) - div( eps0*E_k+1 ) - div( P_k+1 )
//     *
//     *      electromagnetics
//     *        res_k+1 = J_lin + J_nonlin( B_k+1 ) - rot( nu0*B_k+1 ) + rot( M_k+1 )
//     *
//     * 2: incremental evaluation
//     *      electrostatics:
//     *        res_k+1 = res_k - f_nonlin( E_k ) + f_nonlin( E_k+1 ) - div( [eps0 + deltaP_k+1/deltaE_k+1]*deltaE_k+1 )
//     *
//     *      electromagnetics:
//     *        res_k+1 = res_k - J_nonlin( B_k ) + J_nonlin( B_k+1 ) - rot( [nu0 + deltaM_k+1/deltaB_k+1]*deltaB_k+1 )
//     */
//    UInt REScase;
//    
//    /*
//     * 0: RHS_k = RES_k
//     * 1: further cases not needed at the moment, but have to be considered, when performing a full stepping
//     *    i.e. we compute the increment towards the solution from the previous timestep instead of the value from the
//     *    last iteration
//     */
//    UInt RHScase;
//    
//    
//    /*
//     * EvalVersion = 0
//     *  > hystoperator does not affect solution process
//     *  > solve linear system with eps0 / nu0
//     *  > activate evaluation of hystoperator during output evaluation
//     *
//     *  for all iterations solve (rhs may have nonlinear, nonhyst parts):
//     *    electrostatics:
//     *
//     *      div( eps0*E^n+1_k ) = f^n+1_k-1
//     *
//     *    electromagnetics
//     *
//     *      rot( nu0*B^n+1_k ) = J^n+1_k-1
//     */
//    if(evalVersion == 0){
//      // DEBUG version; hysteresis not considered during solution process
//      tensorToReturnLHS = 2; //always solve with eps0/nu0
//      tensorToReturnRHS = 2; //eps0/nu0
//      
//      tensorToAddLHS = 1; //does not matter here
//      tensorToAddRHS = tensorToAddLHS; //does not matter here but should be same as for LHS for reassemble flag to become false
//      useLastTS = false; //does not matter here
//      
//      vectorToReturn = 2; //hystOperator will return 0 instead of P/M;
//      REScase = 1; //incremental computation of res makes no sense as we have no deltaMat at all
//      RHScase = 0; //RHS = RES
//      
//    }
//    /*
//     * EvalVersion = 1
//     *  > during first iteration solve:
//     *    electrostatics:
//     *
//     *      div( eps0*deltaE_1 ) = f^n+1_0 - div( eps0*E^n+1_0 ) - div( P^n+1_0 )
//     *
//     *        deltaE_1 = E^n+1_1 - E^n+1_0
//     *
//     *    electromagnetics
//     *
//     *      rot( nu0*deltaB_1 ) = J^n+1_0 - rot( nu0*B^n+1_0 ) + rot( M^n+1_0 )
//     *
//     *        deltaB_1 = B^n+1_1 - B^n+1_0
//     *
//     *  > during later iterations solve:
//     *    electrostatics:
//     *
//     *      div( [eps0 + deltaP_k-1/deltaE_k-1]*deltaE_k ) = f^n+1_k-1 - div( eps0*E^n+1_k-1 ) - div( P^n+1_k-1 )
//     *
//     *        deltaE_k = E^n+1_k - E^n+1_k-1
//     *
//     *    electromagnetics
//     *
//     *      rot( [nu0 + deltaM_k-1/deltaB_k-1]*deltaB^k ) = J^n+1_k-1 - rot( nu0*B^n+1_k-1 ) + rot( M^n+1_k-1 )
//     *
//     *      deltaB_k = B^n+1_k - B^n+1_k-1
//     *
//     */
//    else if(evalVersion == 1){
//      // DELTA_MAT; rhs and lhs are affected by output of hysteresis operator
//      //  during first iteration use no deltaMatrix
//      if(callingCnt == 0){
//        tensorToReturnLHS = 5; // 2 for eps0/nu0 // 1 to start with epsInitial instead // 4 start with estimated eps
//        tensorToAddLHS = 2; // not needed here
//      } else {
//        tensorToReturnLHS = 3; // deltaMatrix
//        tensorToAddLHS = 2; //add eps0 / nu0
//      }
//      
//      useLastTS = false;
//      if(forceReevaluation || (callingCnt == 0)){
//        // full evaluation of residual; matrix has to be assemble with eps0/nu0; P/M has to
//        // be put on rhs
//        tensorToReturnRHS = 2; //eps0
//        vectorToReturn = 1; //return current state of hystOperator (P/M)
//        tensorToAddRHS = 1; //does not matter here
//        REScase = 1;
//      } else {
//        tensorToReturnRHS = 3; //assemble deltaMat
//        vectorToReturn = 2; //P is already included in deltaMat > return 0 instead
//        tensorToAddRHS = 2; //assemble matrix with eps0 added to deltaP/deltaE / nu0 added deltaM/deltaB
//        REScase = 2; //incremental residual
//      }
//      
//      RHScase = 0; // RHS = RES
//    }
//    /*
//     * EvalVersion = 2
//     *  > during all iterations solve:
//     *    electrostatics:
//     *
//     *      div( [eps0 + deltaP_k-1/deltaE_k-1]*deltaE_k ) = f^n+1_k-1 - div( eps0*E^n+1_k-1 ) - div( P^n+1_k-1 )
//     *
//     *        deltaE_k = E^n+1_k - E^n+1_k-1
//     *        deltaE_0 = E^n+1_0 - E^n  > useLastTS = true
//     *
//     *    electromagnetics
//     *
//     *      rot( [nu0 + deltaM_k-1/deltaB_k-1]*deltaB^k ) = J^n+1_k-1 - rot( nu0*B^n+1_k-1 ) + rot( M^n+1_k-1 )
//     *
//     *      deltaB_k = B^n+1_k - B^n+1_k-1
//     *      deltaB_0 = B^n+1_0 - B^n    > useLastTS = true
//     *
//     */
//    else if(evalVersion == 2){
//      // DELTA_MAT; rhs and lhs are affected by output of hysteresis operator
//      tensorToReturnLHS = 3; // deltaMatrix
//      tensorToAddLHS = 2; //add eps0 / nu0
//      
//      if(callingCnt == 0){
//        useLastTS = true;
//      } else {
//        useLastTS = false;
//      }
//      
//      if(forceReevaluation || (callingCnt == 0)){
//        // full evaluation of residual; matrix has to be assemble with eps0/nu0; P/M has to
//        // be put on rhs
//        tensorToReturnRHS = 2; //eps0
//        vectorToReturn = 1; //return current state of hystOperator (P/M)
//        tensorToAddRHS = 1; //does not matter here
//        REScase = 1;
//      } else {
//        tensorToReturnRHS = 3; //assemble deltaMat
//        vectorToReturn = 2; //P is already included in deltaMat > return 0 instead
//        tensorToAddRHS = 2; //assemble matrix with eps0 added to deltaP/deltaE / nu0 added deltaM/deltaB
//        REScase = 2; //incremental residual
//      }
//      
//      RHScase = 0; // RHS = RES
//    }
//    else{
//      EXCEPTION("Evalversion > 3 not implemented yet");
//    }
//    
//    if((tensorToReturnRHS == 3)&&(tensorToReturnLHS == 3)){
//      /*
//       * in case of delta mat on lhs and rhs, check also if the addTensor variable is the same
//       */
//      reassemble = (tensorToAddLHS != tensorToAddRHS);
//    } else {
//      reassemble = (tensorToReturnLHS != tensorToReturnRHS);
//    }
//    
//    /*
//     * #### PREPARATIONS ####
//     */
//    /*
//     * set current solution by overwriting solVec_
//     */
//		
//    //		std::cout << "oldSolution: " << oldSolution.ToString() << std::endl;
//    //		std::cout << "solVec_.ToString: "  << solVec_.ToString() << std::endl;
//		
//    if(usedEta != 0){
//      //      std::cout << "solVec before add " << solVec_.ToString() << std::endl;
//      //      std::cout << "oldSolution before add " << oldSolution.ToString() << std::endl;
//      //      std::cout << "solIncrement before add " << solIncrement.ToString() << std::endl;
//      //      std::cout << "usedEta " << usedEta << std::endl;
//      
//      solVec_.Add( 1.0, oldSolution, usedEta, solIncrement);
//      //      std::cout << "solVec after add " << solVec_.ToString() << std::endl;
//      //      std::cout << "oldSolution after add " << oldSolution.ToString() << std::endl;
//      //      std::cout << "solIncrement after add " << solIncrement.ToString() << std::endl;
//      //      std::cout << "usedEta " << usedEta << std::endl;
//      
//    } else {
//      solVec_ = oldSolution;
//    }
//    
//    /*
//     * we have to set the stageSol to the actual solutionVector;
//     * this is needed for the handling of the time stepping scheme
//     * (see further below)
//     * by doing so, we can tell the ComputeStageRHS function to include
//     * the current solution and thus consider an incremental stepping
//     * 
//     * Attention: this is only valid if we really want to subtract the whole
//     * solution from the rhs! incremental computation of residual is not possible
//     * thus!
//     */
//    stageSol = solVec_;
//		
//    SBM_Vector outputRHS(BaseMatrix::DOUBLE);
//    if(debugOutput){
//      algsys_->GetRHSVal( outputRHS );
//      //      std::cout << "RHS before setting up residual: " << outputRHS.ToString() << std::endl;
//      //
//      //      std::cout << "old solution (oldSol): " << oldSolution.ToString() << std::endl;
//      //      std::cout << "solInc (solInc): " << solIncrement.ToString() << std::endl;
//      //      std::cout << "usedEta: " << usedEta << std::endl;
//      //      std::cout << "current solution (oldSol + eta*solInc): " << solVec_.ToString() << std::endl;
//      //
//      //      std::cout << "Input parameter:" << std::endl;
//      //      std::cout << "evalVersion: " << evalVersion << std::endl;
//      //      std::cout << "forceReevaluation: " << forceReevaluation << std::endl;
//      //      std::cout << "callingCnt: " << callingCnt << std::endl;
//      //
//      //      std::cout << "Derived flags:" << std::endl;
//      //      std::cout << "REScase: " << REScase << std::endl;
//      //
//      //      std::cout << "tensorToReturnLHS: " << tensorToReturnLHS << std::endl;
//      //      std::cout << "tensorToReturnRHS: " << tensorToReturnRHS << std::endl;
//      //
//      //      std::cout << "tensorToAddLHS: " << tensorToAddLHS << std::endl;
//      //      std::cout << "tensorToAddRHS: " << tensorToAddRHS << std::endl;
//      //
//      //      std::cout << "vectorToReturn: " << vectorToReturn << std::endl;
//      //      std::cout << "useLastTS: " << useLastTS << std::endl;
//      //
//      //      std::cout << "Reassembling needed: " << reassemble << std::endl;
//      
//    }
//    
//    /*
//     * evaluate and store current linear RHS
//     * > has to be done only once for each TS > during first call
//     */
//    if(callingCnt == 0){
//      if(debugOutput){
//        std::cout << "Assemble lin rhs" << std::endl;
//      }
//      currentLinRhsVec_.Init();
//      algsys_->InitRHS(currentLinRhsVec_); //load storage into algsys
//      assemble_->AssembleLinRHS(); //write current linear rhs to algsys
//      algsys_->GetRHSVal( currentLinRhsVec_ ); //retrieve value
//    }
//    
//    //				std::cout << "currentLinRhsVec_ " << currentLinRhsVec_.ToString() << std::endl;
//    //		
//		
//    /*
//     * evaluate and store current NON-linear RHS (without hysteresis!)
//     * > has to be done during each iteration as it depends on the current solution
//     */
//    if(debugOutput){
//      std::cout << "Assemble nonlin rhs (without hysteresis)" << std::endl;
//    }
//    currentNonLinRhsVec_.Init();
//    algsys_->InitRHS(currentNonLinRhsVec_); //load storage into algsys
//    PDE_.SetFlagInCoefFncHyst("vectorToReturn",2); // deactivate hysteresis!
//    assemble_->AssembleNonLinRHS(); //write current linear rhs to algsys
//    algsys_->GetRHSVal( currentNonLinRhsVec_ ); //retrieve value
//    
//    //		std::cout << "currentNonLinRhsVec_ " << currentNonLinRhsVec_.ToString() << std::endl;
//    //		
//		
//    /*
//     * #### EVALUATE RESIDUAL ACCORDING TO PREVIOUSLY SET FLAGS ####
//     */
//    /*
//     * Assemble matrices
//     */
//    // setup flags for matrix computation
//    PDE_.SetFlagInCoefFncHyst("tensorToReturn",tensorToReturnRHS);
//    PDE_.SetFlagInCoefFncHyst("tensorToAdd",tensorToAddRHS);
//    PDE_.SetFlagInCoefFncHyst("useLastTS",useLastTS);
//    assemble_->AssembleMatrices(false);
//    
//    /*
//     * Setup residual
//     */
//    currentResVec_.Init();
//    
//    /* REScase:
//     * 1: full evaluation
//     *        res_k+1 = f_lin + f_nonlin( E_k+1 ) - div( eps0*E_k+1 ) - div( P(E_k+1) )
//     *
//     * 2: incremental computation
//     *        res_k+1 = res_k + f_nonlin( E_k+1 ) - f_nonlin( E_k ) - div( deltaMat_k*deltaE_k )
//     *
//     */
//    debugOutput = true;
//    if(REScase == 1){
//      if(debugOutput){
//        std::cout << "RES - full reevaluation" << std::endl;
//      }
//      
//      /* REScase:
//       * 1: full evaluation
//       *        res_k+1 = f_lin + f_nonlin( E_k+1 ) - div( P(E_k+1) ) - div( eps0*E_k+1 )
//       */
//      // add f_lin
//      currentResVec_.Add(1.0, currentLinRhsVec_);
//      
//      // put res to rhs
//      algsys_->InitRHS(currentResVec_);
//      
//      if(debugOutput){
//        algsys_->GetRHSVal( outputRHS );
//        std::cout << "RES - linear part: " << outputRHS.ToString() << std::endl;
//      }
//      
//      // TODO clean up!
//      // temporary fix
//      // for full reevaluation, we have to evaluate the rhs loading with
//      // the current state of the hyst operator
//      PDE_.SetFlagInCoefFncHyst("forceCurrentTS",1); 
//
//      
//      // + f_nonlin( E_k+1 ) - div( P(E_k+1) )
//      PDE_.SetFlagInCoefFncHyst("vectorToReturn",vectorToReturn); // activate hysteresis!
//      assemble_->AssembleNonLinRHS(); //add current nonlinear rhs to algsys
//      
//      // reset
//      PDE_.SetFlagInCoefFncHyst("forceCurrentTS",0); 
//      
//      if(debugOutput){
//        algsys_->GetRHSVal( outputRHS );
//        std::cout << "RES - linear+nonlinear+hyst part: " << outputRHS.ToString() << std::endl;
//      }
//			
//      //			std::cout << "solVec_ "  << solVec_.ToString() << std::endl;
//      //			std::cout << "solVec_[0].GetSize() " << solVec_(0).GetSize() << std::endl;
//      //      algsys_->GetRHSVal( outputRHS );
//      //			std::cout << "outputRHS "  << outputRHS.ToString() << std::endl;
//      //			std::cout << "outputRHS[0].GetSize()" << outputRHS(0).GetSize() << std::endl;
//      
//      // - div( eps0*E_k+1 )
//      solVec_.ScalarMult(-1.0);
//			
//			algsys_->UpdateRHS(STIFFNESS,solVec_,true);
//			
//			if(!trans){
//				algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
//			}
//      //			
//      //			if(!trans){
//      //				algsys_->UpdateRHS(STIFFNESS,solVec_,true);
//      //			} else {
//      //				// for systems without mass/damping etc it was enough to just update
//      //				// stiffness and stiffness update
//      //				// howerver; for piezo, magnetic and magstrict, we have to update with
//      //				// mass, damping etc, too
//      //				// otherwise, the contributions from those matrices won't get to 0
//      //			 // if(trans){
//      //			 //   algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
//      //			 // }
//      //				// i.e. we have to update rhs with Keff*oldSolution or Keff*currentSolution
//      //				/*
//      //				 * However, the following lines do not work as intended, as the do not
//      //				 *					compute Keff but Kmod = K*D*M instead
//      //				 *				 (i.e. the time integration factors are missing in front of
//      //				 *					D and M)
//      //				 * However 2, this proceducre might not be needed at all as the
//      //				 *				 there is the possibility to add the current solution
//      //				 *				 directly to the computation of the stage RHS
//      //				 *				 see TimeSchemeGLM.cc   void TimeSchemeGLM::ComputeStageRHS(UInt actStage, Integer derivId,
//      //                                      SingleVector* rhsVec, Integer subIdx)
//      //				 *				According to this function, the stage RHS will be considered
//      //				 *				automatically if the nonlin stepping type is INCREMENTAL
//      //				 * Open Question:
//      //				 *				Does it only update with damping and mass or also with stiffness?
//      //				 *
//      //				 */
//      //				
//      //				std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
//      //				std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
//      //				std::map<FEMatrixType,Integer>::iterator matIt;
//      //	//      UInt pos = 0;
//      //
//      //				// update RHS with all matrices!
//      //				for(matIt = matrices.begin();matIt != matrices.end();matIt++){
//      //					if(matIt->second < 0)
//      //						continue;
//      //
//      //					algsys_->UpdateRHS(matIt->first,solVec_,true);
//      //				}
//      //			}
//      
//      solVec_.ScalarMult(-1.0);
//      
//      if(debugOutput){
//        algsys_->GetRHSVal( outputRHS );
//        std::cout << "RES - linear part+nonlinear+hyst part - solution: " << outputRHS.ToString() << std::endl;
//        std::cout << "solution: " << solVec_.ToString() << std::endl;
//        std::cout << "solIncrement: " << solIncrement.ToString() << std::endl;
//      }
//    } else if (REScase == 2){
//      EXCEPTION("Incremental computation not fixed yet");
////			
////			if(debugOutput){
////        std::cout << "RES - incremental evaluation" << std::endl;
////        std::cout << "oldItResVec: " << oldItResVec_.ToString() << std::endl;
////      }
////      
////      /* REScase:
////       * 4: incremental computation
////       *        res_k+1 = res_k + f_nonlin( E_k+1 ) - f_nonlin( E_k ) - div( deltaMat_k*deltaE_k )
////       */
////      // add res_k
////      currentResVec_.Add(1.0, oldItResVec_);
////      
////      // + f_nonlin( E_k+1 ) - f_nonlin( E_k )
////      currentResVec_.Add(1.0, currentNonLinRhsVec_);
////      currentResVec_.Add(-1.0, oldItNonLinRhsVec_);
////      
////      // put res to rhs of system
////      algsys_->InitRHS(currentResVec_);
////      
////      //- div( deltaMat_k*deltaE_k )
////      solIncrement.ScalarMult(-1.0*usedEta);
////			
////			// test: update rhs with keff*solincrement
////			// however this might not be possible as we do not get the 
////			// correct scalings for damping and mass matrix 
////			// i.e. by the following, we apply
////			// (K+M*D)*solIncrement instead of
////			// (K+1/beta*dt^2 M + gamma/dt D)*solIncrement
////			// > here we might need an additional function that
////			// sclaes the solution vector by the correct entries
////			
////			
////			
////			if(!trans){
////				algsys_->UpdateRHS(STIFFNESS,solIncrement,true);
////			} else {
////				// for systems without mass/damping etc it was enough to just update
////				// stiffness and stiffness update
////				// howerver; for piezo, magnetic and magstrict, we have to update with
////				// mass, damping etc, too
////				// otherwise, the contributions from those matrices won't get to 0
////        // if(trans){
////        //   algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
////        // }
////				// i.e. we have to update rhs with Keff*oldSolution or Keff*currentSolution
////				std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
////				std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
////				std::map<FEMatrixType,Integer>::iterator matIt;
////        //      UInt pos = 0;
////        
////				// update RHS with all matrices!
////				for(matIt = matrices.begin();matIt != matrices.end();matIt++){
////					if(matIt->second < 0)
////						continue;
////          
////					algsys_->UpdateRHS(matIt->first,solIncrement,true);
////				}
////			}
////			
////      //			// old:
////      //      algsys_->UpdateRHS(STIFFNESS,solIncrement,true);
////      //
////      //      if(trans){
////      //        algsys_->UpdateRHS(STIFFNESS_UPDATE,solIncrement,true);
////      //      }
////      
////      solIncrement.ScalarMult(-1.0/usedEta);
////      
//    } else {
//      EXCEPTION("REScase unknown");
//    }
//    
//    //retrieve value
//    algsys_->GetRHSVal( currentResVec_ );
//    
//    if(debugOutput){
//      algsys_->GetRHSVal( outputRHS );
//      std::cout << "RES after setting up residual; before considering timestepping: " << outputRHS.ToString() << std::endl;
//    }
//    
//    /*
//     * common steps and all methods and callingCnt
//     */
//    /*
//     * update with stage rhs / updating accoring to time stepping
//     */
//    if(trans){
//      
//      std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
//      std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
//      std::map<FEMatrixType,Integer>::iterator matIt;
//      //      UInt pos = 0;
//      
//      // update RHS according to time stepping
//			// in order to correctly subtract the current solution from the rhs
//			// we have to make sure that the nonLinType is INCREMENTAL
//			// If this is the case, ComputeStageRHS will also add the correctly
//			// scaled stage solution to the stageRHS and thus we get
//			// -(Keff-K)*stageSol on the rhs (note: K is usually not included in ComputeStageRHS
//			// as the corresponding factor = 0 (for Newmark); so we have to subtract
//			// K*stageSol separately (see above)
//			bool forceIncremental = true;
//      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
//        if(matIt->second < 0)
//          continue;
//        
//        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
//          FeFctIdType fctId = fncIt->second->GetFctId();
//          fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(fctId),-1,forceIncremental);
//        }
//        algsys_->UpdateRHS(matIt->first,stageRHS_,true);
//      }
//			
//			// same issue maybe as with the initialisation > for coupled pde, pos does not give
//			// the right subindex
//      //      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
//      //
//      //        if(matIt->second < 0)
//      //          continue;
//      //        for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//      //          fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(pos));
//      //        }
//      //        algsys_->UpdateRHS(matIt->first,stageRHS_,true);
//      //      }
//    }
//    
//    if(debugOutput){
//      algsys_->GetRHSVal( outputRHS );
//      std::cout << "RES after setting up residual; after considering timestepping: " << outputRHS.ToString() << std::endl;
//    }
//    debugOutput = false;
//    /*
//     * update current residual with information of time stepping
//     * TODO: not actually sure about this step; should the time stepping update count to residual or not?
//     *       pro:
//     *        via this step we take mass matrix * second time derivative into account
//     *        including these values will actually lead to decreasing residual which is not the case without the addition
//     *
//     */
//    algsys_->GetRHSVal( currentResVec_ );
//    
//        std::cout << "reassemble: " << reassemble << std::endl;
//        std::cout << "skipEeasembly: " << skipReassembly << std::endl;
//    /*
//     * #### REASSEMBLE SYSTEM ####
//     */
//    if((reassemble)&&(skipReassembly == false)){
//      
//      if(debugOutput){
//        std::cout << "Re-Assembling system for next iteration" << std::endl;
//      }
//      
//      //        std::cout << "tensorToReturnLHS: " << tensorToReturnLHS << std::endl;
//      PDE_.SetFlagInCoefFncHyst("tensorToReturn",tensorToReturnLHS);
//      PDE_.SetFlagInCoefFncHyst("tensorToAdd",tensorToAddLHS);
//      PDE_.SetFlagInCoefFncHyst("useLastTS",useLastTS);
//      assemble_->AssembleMatrices(false);
//    }
//    					
//				
//    /*
//     * #### REPLACE RHS (full stepping only) ####
//     */
//    if(RHScase != 0){
//      EXCEPTION("Currently only incremental stepping implemented");
//    }
//    
//  }
//    
//  Double StdSolveStep::LineSearchHyst(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double& etaLineSearch, UInt evalVersion,
//          UInt callingCnt, bool trans, bool performLineSearch,
//          bool forceReevaluation, bool debugOutput, bool reset, UInt allowedSteps)  {
//    
//    static Timer* lineSearchTimer = new Timer();
//    static UInt numCalls = 0;
//    if(debugOutput){
//      std::cout << "LineSearchHyst" << std::endl;
//      lineSearchTimer->Start();
//      numCalls++;
//    }
//    
//    /*
//     * backup nonlinear rhs, residual and solution vector
//     */
//    SBM_Vector solOld(BaseMatrix::DOUBLE);
//    solOld = solVec_;
//    
//    SBM_Vector resOld(BaseMatrix::DOUBLE);
//    resOld = currentResVec_;
//    
//    
//    const UInt nrEtas = 13;
//    //Double eta[nrEtas] = {1.0,0.5,0.2,0.1};
//    Double eta[nrEtas] = {1.0,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.0025,0.001,0.0005,0.00025,0.0001};
//    //Double eta[nrEtas] = {1.0,0.316,0.1,0.0316,0.01,0.00316,0.001,0.000316,0.0001,0.0000316,0.00001,0.0000316,0.00001};
//    //Double eta[nrEtas] = {1.0,0.075,0.05,0.025,0.1};
//    //Double eta[nrEtas] = {1.0,0.316,0.1,0.0316,0.01};
//    //Double eta[nrEtas] = {1.0,0.75,0.5,0.25,0.1};
//    bool assumeQuadratic = false;
//    
//    // initialize etaOpt or receive compiler warning
//    Double etaOpt = 0.0;
//    Double residualL2NormOpt = 1e15;
//    
//    if(performLineSearch == false){
//      /*
//       * only take first eta = 1
//       * -> we do not have to lock hysteresis nor do we have to backup
//       *    the solution or the residual as we directly take the first value
//       */
//      etaOpt = 1.0;
//    } else {
//      /*
//       * lock hysteresis memory
//       * NOTE: in the current evaluation, the locking is done automatically
//       * in coefFnc hyst if the system is called during assembly;
//       * only if SetPreviousHystValues is used with setLastTS = true, the
//       * hysteresis memory is unlocked
//       */
//      
//      /*
//       * test different eta
//       */
//      bool skipReassembly = true;
//      UInt nrEtasToTest = std::min(allowedSteps,nrEtas);
//      for( UInt i=0; i<nrEtasToTest; i++) {
//        if(debugOutput){
//          std::cout << "testing eta: " << eta[i] << std::endl;
//        }
//        
//        CalcResidualAndConfigSystemForHysteresis(solOld,solIncrement,stageSol,eta[i], 0, callingCnt, evalVersion,
//                trans, forceReevaluation, skipReassembly, debugOutput, reset);
//        
//        /*
//         * calculate error norm as criterion
//         */
//        Double residualL2Norm = currentResVec_.NormL2();
//        
//        if(debugOutput){
//          std::cout << "L2-norm of computed residual for eta = " << eta[i] << ": " << residualL2Norm << std::endl;
//        }
//        
//        /*
//         * restore old solution vector, old residual and nonlinear rhs
//         * -> needed to perform a correct next call to CalcResidual...
//         */
//        solVec_ = solOld;
//        currentResVec_ = resOld;
//        
//        if (residualL2Norm < residualL2NormOpt) {
//          residualL2NormOpt = residualL2Norm;
//          etaOpt = eta[i];
//        } else {
//          if(assumeQuadratic){
//            /*
//             * here we assume that the residual follows a quadratic curve
//             * i.e. we only have one minimum which is the global minimum;
//             * in that case, if the residual increases between two iterations,
//             * we can assume that we passed the minimum already and stop
//             */
//            break;
//          }
//        }
//      }
//      
//      /*
//       * unlock hysteresis memory
//       * -> no longer done here but first when we store the lastTS values
//       */
//      
//    }
//    
//    /*
//     * now we have found etaOpt, so we can call CalcResidualAndConfigSystemForHysteresis
//     * one last time to setup the correct system for next computation
//     */
//    etaLineSearch = etaOpt;
//    
//    if(debugOutput){
//      lineSearchTimer->Stop();
//      std::cout << "Total time needed for linesearch: " << lineSearchTimer->GetCPUTime() << std::endl;
//      std::cout << "Average time needed for linesearch: " << lineSearchTimer->GetCPUTime()/numCalls << std::endl;
//      std::cout << "(optimal) eta: " << etaOpt << std::endl;
//    }
//    
//    // calculate res and assemble system once more using the optimal eta
//    bool skipReassembly = false;
//    CalcResidualAndConfigSystemForHysteresis(solOld,solIncrement,stageSol,etaLineSearch, 0, callingCnt, evalVersion,
//            trans, forceReevaluation, skipReassembly, debugOutput, reset);
//    
//    if(debugOutput){
//      std::cout << "minimal residual: " << residualL2NormOpt << std::endl;
//    }
//    
//    return currentResVec_.NormL2();
//  }
//  
//  void StdSolveStep::StepTransNonLinHysteresis() {
//    /*
//     * New approach / new hope?
//     *
//     * Idea:
//     *  1. Iteration: full step starting at 0
//     *    i.e. Reset solution vector to 0
//     *         Reset old iteration vector to 0
//     *         Lock memory of hysteresis operator
//     *          > it will act as a nonlinear curve; no branching
//     *          > due to reset of solution to 0 it will be dropped from its last state
//     *            (= solution state after last timestep) to a remanent state
//     *            > this drop is not stored in the memory of the operator, i.e.
//     *              if we would apply the solution of the last ts to the operator, we
//     *              would end up at the same state as prior to the drop to remanence
//     *         Assemble system using the difference between the values from lastTS
//     *         and the reset solution
//     *  2-n. Iterartion: incremental stepping
//     *    i.e. Still lock hyst operator
//     *         Assemble system using difference from current state to the lastIT state
//     * */
//    static UInt timestepCnt = 0;
//    timestepCnt++;
//    if(debugOutput_){
//      std::cout << "####                   #### " << std::endl;
//      std::cout << "#### BEGIN OF TIMESTEP #### " << std::endl;
//      std::cout << "####   Nr: " << timestepCnt << "   #### " << std::endl;
//    }
//    
//    LOG_TRACE(stdsolvestep) << "StdSolveStep::StepTransNonLinHysteresis";
//    
//    SBM_Vector solInc(BaseMatrix::DOUBLE);
//    SBM_Vector actSol(BaseMatrix::DOUBLE);
//    
//    //		std::cout << "solVec_: " << solVec_.ToString() << std::endl;
//    //		
//    bool deltaIDBC = true;
//    actSol = solVec_;
//		
//    //obtain the number of stages
//    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
//    if ( numStages > 1 ){
//      /*
//       * maybe it is quite easy to enable multiple stages, maybe it does not work;
//       * feel free to find out
//       */
//      EXCEPTION("StepTransNonLinHysteresis: just one stage time-stepping allowed");
//    }
//    
//    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
//    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
//    std::map<FEMatrixType,Integer>::iterator matIt;
//    
//    UInt pos = 0;
//    
//    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
//    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
//      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
//    }
//    
//    //do initialization
//    rhsVec_.Init();
//    algsys_->InitRHS(rhsVec_);
//    
//    LOG_DBG(stdsolvestep) << "StepTransNonLinHysteresis: Stage: " << 0 ;
//    
//    //we obtain a reference to the stage vectors of the scheme
//    /*
//     * Important: Although stageSol seems not to be necessary, the TimeScheme will rely on it
//     * It HAS to store the total solution (not the increment!) -> see TimeSchemeGLM.cc
//     * Important2: Incremental stepping is only working (correctly) for EffectiveStiffness formulation
//     */
//    SBM_Vector stageSol;
//    stageSol.Resize(feFunctions_.size());
//		
//    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
//      FeFctIdType fctId = fncIt->second->GetFctId();
//      stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(0),fctId);
//      fncIt->second->GetTimeScheme()->InitStage(0,actTime_,PDE_.GetDomain());
//    }
//		
//		// so nicht! dadurch ist die stageSol womöglich falsch aufgebaut (z.b. mech und
//		// elec unknows vertauscht!
//    //    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//    //      stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(0),pos);
//    //      fncIt->second->GetTimeScheme()->InitStage(0,actTime_,PDE_.GetDomain());
//    //    }
//    stageSol.SetOwnership(false);
//		
//    /*
//     * During the later computations, we will add all incremental solutions to
//     * stageSol; unlike the case of StepTransNonLin, we add the results to the
//     * current solution instead of only gathering the increments.
//     * Note that we initialize stageSol with actSol in StepTransNonLin, too, but
//     * there, we will Init stageSol with 0.0 during the first iteration
//     */
//    stageSol = actSol;
//    solVec_  = actSol;
//    
//		/*
//     Irgendetwas stimmt hier doch nicht!
//     Am Anfang hat solVec_ die richtige Größe
//     Nachdem actSol an stageSol zugewiesen wurde, hat solVec plötzlich 
//     subVector 1 und 2 vertauscht, stageSol und actSol sind aber noch richtig
//     
//     Danach nachdem solVec = actSol passen die Größen von solVec und actSol
//     wieder, aber stageSol hat jetzt vertauschte einträge
//     * 
//     * Setzt man statt solVec = actSol, solVec = stageSol so passt gar nichts
//     * mehr
//     */
//		
//    // setup right hand side
//    Double loadFactor = 1.0;
//    Double RhsLinL2Norm = SetLinRHS(loadFactor);
//    
//    /*
//     * disable output of switching and rotation state as BMP images (even if flag
//     * printOut is set to 1 in material file)
//     */
//    PDE_.SetFlagInCoefFncHyst("allowBMP",false);
//    
//    /*
//     * START ITERATIONS
//     */
//    static UInt overallIterationCounter = 0;
//    bool performOneMoreStep;
//    bool transient = true;
//    UInt iterationCounter = 0;
//    UInt stage = 0;
//    Double initialEta = 0.0;
//    Double residualErr;
//    Double incrementalErr;
//    bool forceReevaluation = forceReevaluation_;
//    
//    /*
//     * NOTE:
//     * evalVersion_
//     * 0 = debug fixpoint > hysteresis will be evaluated only during output
//     *
//     * 1 = fixpoint / deltaMat
//     *       lhs and rhs with eps0/nu0 > only during first iteration!
//     *       lhs = deltaMat, rhs with eps0/nu0 during later iteration
//     * 2 = deltaMat during all iterations
//     */
//    //    std::cout << "evalVerion: " << evalVersion_ << std::endl;
//    
//    /*
//     * ITERATION 0
//     * > Prepare initial system
//     */
//    PDE_.SetFlagInCoefFncHyst("estimateSlope",9);
//    bool skipReassembly = false; // currently not used
//    
//    bool reset = false; 
//    UInt numResets = 0;
//    
//    // keep track of residual; if it is not decreasing over several iterations > reset (see below)
//    std::deque<Double> resNormHistory(5, 50000);
//    
//    std::cout << "solVec_ at start of TS: " << solVec_.ToString() << std::endl;
//    std::cout << "stageSol at start of TS: " << stageSol.ToString() << std::endl;
//    std::cout << "actSol at start of TS: " << actSol.ToString() << std::endl;
//    
//    CalcResidualAndConfigSystemForHysteresis(solVec_,solInc,stageSol, initialEta, stage,
//            iterationCounter, evalVersion_, transient,
//            forceReevaluation, skipReassembly, debugOutput_, reset);
//    // after this initial setup, we have to call SetLastItOrLastTSSBMVectors, to store the current
//    // res and nonlinRHS as the values from the previous iteration so that these values are correct
//    // when we call CalcResidualAndConfigSystemForHysteresis during the first actual iteration
//    SetLastItOrLastTSSBMVectors(false);
//    /*
//     * ITERATION > 0
//     */
//    do {
//      
//      overallIterationCounter++;
//      iterationCounter++;
//      
//      LOG_DBG(stdsolvestep) << "StepTransNonLinHysteresis: Iteration: " << iterationCounter ;
//      
//      if(debugOutput_){
//        std::cout << "####                   #### " << std::endl;
//        std::cout << "#### BEGIN OF ITERATION #### " << std::endl;
//        std::cout << "####   Nr: " << iterationCounter << "   #### " << std::endl;
//        std::cout << "####   Total Nr: " << overallIterationCounter << "   #### " << std::endl;
//      }
//      
//      /*
//       *  Assembly has already be done during CalcResidualAnd...
//       *  Now build up the algebraic system
//       */
//      // set system matrix to zero initially, as ConstructEffectiveMatrix only
//      // sums up the contributions
//      matrix_factor_.clear();
//      algsys_->InitMatrix(SYSTEM);
//      for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
//        FeFctIdType fctId = fncIt->second->GetFctId();
//        fncIt->second->GetTimeScheme()->AddMatFactors(0,matrices,matrix_factor_[fctId]);
//        algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
//      }
//      
//      PDE_.SetBCs();
//      algsys_->BuildInDirichlet();
//      algsys_->SetupPrecond();
//      algsys_->SetupSolver();
//      
//      bool setIDBC = false;
//      
//      if (( iterationCounter == 1 && couplingIter_ == 0 )||(reset)) {
//        /*
//         * we only need to include the idbc values once during the first iteration
//         * in the new approach, we do not need deltaIDBC as we start from 0
//         */
//        setIDBC = true;
//      }
//      
//      /*
//       *  SOLVE SYSTEM
//       */
//      algsys_->Solve(setIDBC,deltaIDBC);
//      
//      /*
//       *  RETRIEVE SOLUTION INCREMENT
//       *  1, iteration: increment towards 0, i.e. full solution
//       *  2+ iteration: increment towards previous iteration
//       */
//      algsys_->GetSolutionVal(solInc, setIDBC, deltaIDBC );
//      
//      /*
//       * CALCULATE RESIDUAL AND SETUP SYSTEM FOR NEXT ITERATION VIA LINESEARCH
//       */
//      Double residualL2Norm = 0.0;
//      Double etaLineSearch  = 1.0;
//      bool performLinesearch;
//      
//      if (( lineSearch_ == "none" || iterationCounter == 1 )||(reset)) {
//        /*
//         * during the first step, we need a full step with eta = 1.0
//         * this has to be done to ensure correct idbc values!
//         */
//        performLinesearch = false;
//      } else {
//        performLinesearch = true;
//      }
//      
//      if(reset){
//        // reset reset flag
//        reset = false;
//      }
//      
//      // the linesearching will test different etas and take the one which
//      // produces the smallest residual
//      // afterwards the system will be setup (assembled) with this new solution
//      // new solution (oldSolution + etaOpt*solutionIncrement) is written to solVec_
//      // with each reset, try also one smaller step for linesearch
//      UInt allowedSteps = 7+2*numResets;
//      residualL2Norm = LineSearchHyst(solInc, stageSol, etaLineSearch, evalVersion_, iterationCounter,
//              transient , performLinesearch ,forceReevaluation, debugOutput_, reset, allowedSteps);
//      
//      // calculation of relative residual error =======================================
//      //if ( RhsLinL2Norm > 0.0 ){
//      //Double oldResError = residualErr;
//      if ( RhsLinL2Norm > 1.0 ){
//        residualErr = residualL2Norm / RhsLinL2Norm;
//      } else {
//        residualErr = residualL2Norm;
//      }
//      
//      resNormHistory.push_back(residualErr);
//      resNormHistory.pop_front();
//      
//      // calculate incremental error ========================================
//      //TODO: we should have eta*solInc.NormL2() as we do not perform full steps in general
//      //Double solIncrL2Norm = solInc.NormL2();
//      
//      Double solIncrL2Norm = etaLineSearch*solInc.NormL2();
//      Double solL2Norm  = solVec_.NormL2();
//      
//      if ( solL2Norm > 1.0 ){
//        // if ( actSolL2Norm > 1e-7 ){
//        incrementalErr = solIncrL2Norm / solL2Norm;
//      } else {
//        incrementalErr = solIncrL2Norm;
//        // incrementalErr = solIncrL2Norm/1e-7;
//      }
//      
//      // update actSol to currently computed value
//      actSol = solVec_;
//      
//      performOneMoreStep =
//              (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
//      
//      if(residualErr < 1e-18){
//        /*
//         * in some cases e.g. E -> 0, P != 0 (remanence)
//         * the residual might become very small although the solution does not converge
//         * if the residual becomes too small ( ~ 1e-28) the system reports
//         * error during solution (at least if using pardiso);
//         * -> a complete reevaluation of the residual might help
//         *
//         */
//        forceReevaluation = true;
//      } else {
//        /*
//         * restore old state
//         */
//        forceReevaluation = forceReevaluation_;
//      }
//      
//      if(debugOutput_){
//        std::cout << "UsedEta: " << etaLineSearch << std::endl;
//        std::cout << "Solution update: " << solInc.ToString() << std::endl;
//        std::cout << "current solution norm: " << solL2Norm << std::endl;
//        std::cout << "incrementalErr_abs: " << solIncrL2Norm << std::endl;
//        std::cout << "incrementalErr_rel: " << incrementalErr << std::endl;
//        std::cout << "residualErr: " << residualErr << std::endl;
//      }
//      
//      // compare residual at front and back; should have decreased; otherwise try reset
//      if(resNormHistory.front() <= resNormHistory.back()){
//        std::cout << "Residual did not decrease. Reset?" << std::endl;
//        std::cout << "Iteration " << iterationCounter << " of timestep " << timestepCnt << std::endl;
//        std::cout << "Number of previous resets: " << numResets << std::endl;
//        reset = true;
//        
//        // reset deque
//        for(UInt i = 0; i < resNormHistory.size(); i++){
//          resNormHistory[i] = 50000;
//        }
//      }
//      
//      if(reset == true){
//        numResets += 1;
//        //      std::cout << "Perform reset" << std::endl;
//        // set solution to 0 i.e. set new starting point for iterations
//        // note: the memory of the hyst operator is not deleted but as we do not
//        //       write it until the end of the timestep, it is basically still the
//        //       state from the start of the timestep
//        solVec_.Init();
//        // estimate slope around this point
//        // estimate using a stepping
//        PDE_.SetFlagInCoefFncHyst("estimateSlope",9+numResets);
//        // for the first iteration after the reset, we have to include the boundary conditions
//        // again and they must not e the deltaBCs but the full boundary conditions (as we start from 0)
//        deltaIDBC = false;
//        
//        // now we have to reassemble the system (we basically repeat all steps till start of loop)
//        // the reset flag will have the same meaning as an iterationCounter of 0 (i.e. full evaluation, computation of rhs, ...)
//        CalcResidualAndConfigSystemForHysteresis(solVec_,solInc,stageSol, initialEta, stage,
//                iterationCounter, evalVersion_, transient,
//                forceReevaluation, skipReassembly, debugOutput_, reset);
//        
//        if(numResets > 4){
//          EXCEPTION("Error cannot be decreases by resetting");
//        }                                                 
//      }
//      
//      
//      // output of norms and data to info.xml
//      if ( nonLinLogging_ == true ) {
//        // get current step
//        UInt actStep = PDE_.GetSolveStep()->GetActStep();
//        
//        if (PDE_.IsIterCoupled()) {
//          WriteNonLinIterToInfoXML(pdename_, couplingIter_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
//        } else {
//          //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
//          WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
//        }
//        // write norm to file
//        logFile_ <<  iterationCounter << "\t"
//                << residualErr << "\t"
//                << incrementalErr << "\t"
//                << etaLineSearch << std::endl;
//      }
//      
//      /*
//       * STORE CURRENT VALUES FOR NEXT ITERATION
//       */
//      bool lastTS = false;
//      // Note: by calling SetPreviousHystVals with flag lastTS = false,
//      //        coefFncHyst will evaluate the hysteresis operator with
//      //        the current state of solVec_; for this evaluation, the memory of
//      //        the hysteresis operator will be locked, i.e. the evaluation will
//      //        lead to reversible changes
//      PDE_.SetPreviousHystVals(lastTS);
//      SetLastItOrLastTSSBMVectors(lastTS);
//      
//    } while(performOneMoreStep && (iterationCounter < nonLinMaxIter_));
//    
//    abortOnMaxIter_ = true;
//    if (performOneMoreStep && (iterationCounter >= nonLinMaxIter_) && abortOnMaxIter_) {
//      EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
//              << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
//              << "' at iteration '" << iterationCounter
//              << "'.\n ==> incremental error: " << incrementalErr
//              << "\n ==> residual error: " << residualErr);
//    }
//    
//    // check solution
//    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator limitFeFctIt;
//    limitFeFctIt = feFunctions_.find(solutionLimit_);
//    if (limitFeFctIt != feFunctions_.end() ) {
//      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
//          //const SingleVector * subv = solVec_.GetPointer(pos);
//          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(pos)));
//          for (UInt j=0; j < dsubVec.GetSize(); j++) {
//            if (dsubVec[j] >= maxValidValue_) {
//              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
//                      "' is larger than the allowed maximum limit set in the XML: "
//                      << maxValidValue_);
//            }
//            if (dsubVec[j] <= minValidValue_) {
//              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
//                      "' is smaller than the allowed minimum limit set in the XML: "
//                      << minValidValue_);
//            }
//          }
//        }
//      }
//    }
//    
//    //update stage
//    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
//      fncIt->second->GetTimeScheme()->FinishStep();
//    }
//    
//    if(debugOutput_){
//      std::cout << "Iterative scheme was successful after " << iterationCounter << " iterations" << std::endl;
//    }
//    
//    /*
//     * STORE CURRENT VALUES FOR NEXT ITERATION
//     */
//    bool lastTS = true;
//    // Note: by calling SetPreviousHystVals with flag lastTS = true,
//    //        coefFncHyst will evaluate the hysteresis operator with
//    //        the final state of solVec_; for this evaluation, the memory of
//    //        the hysteresis operator will be unlocked, i.e. the evaluation will
//    //        lead to irreversible changes
//    PDE_.SetPreviousHystVals(lastTS);
//    SetLastItOrLastTSSBMVectors(lastTS);
//    
//    /*
//     * Store IDBC values from the current time step in form of its rhs representation
//     * (i.e. the effect that it will have on the rhs)
//     * these values are needed to compute the deltaIDBC values
//     * (currently not used)
//     */
//    algsys_->SetOldDirichletValues();
//    
//    /*
//     * allow to output switching and rotation state as BMP images (if flag
//     * printOut is set to 1 in material file)
//     */
//    PDE_.SetFlagInCoefFncHyst("allowBMP",true);
//    
//    /*
//     * set evaluationPurpose to 4, i.e. output
//     * this will lock the hysteresis operator again and will only evaluate
//     * the hystOperator at the midpoint of each element regardless of the
//     * evaluation depth
//     * > for more infos see coefFncHyst
//     */
//    PDE_.SetFlagInCoefFncHyst("evaluationPurpose",4);
//    
//    if(debugOutput_){
//      std::cout << "####                 #### " << std::endl;
//      std::cout << "#### END OF TIMESTEP #### " << std::endl;
//      std::cout << "####                 #### " << std::endl;
//      PDE_.SetFlagInCoefFncHyst("outputDebugInfos",1);
//    }
//    
//    //TODO: continue here
//    //TODO: add defaultValue to xmlSchema!!!!!
//  }
 
  
  void StdSolveStep::StepTransNonLinMaterial() {
    
    REFACTOR;
    //
    //    bool performOneMoreStep;
    //
    //    SBM_Vector solInc, actSol;
    //    actSol.Init();
    //
    //    // set iteration counter
    //    UInt iterationCounter=0;
    //    Double RhsLinL2Norm;
    //    SBM_Vector uOld, actRHS;
    //
    //    StepTransLin;
    //    algsys_->GetSolutionVal( actSol );
    //    PDE_.SaveSolution( actSol );
    //
    //    // to incorporate loads
    //    Double loadFactor = 1.0;
    //    RhsLinL2Norm = SetLinRHS(loadFactor);
    //
    //
    //    do {
    //      uOld=actSol;
    //      // compute u_{n+1}^k+1
    //      iterationCounter++;
    //
    //      PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(analysis_id, "nonLin", iterationCounter);
    //
    //      // re initialize RHS and system matrix
    //      algsys_->InitRHS();
    //
    //      assemble_->AssembleLinRHS();
    //
    //      assemble_->AssembleMatrices();
    //
    //      // account for Dirichlet BCs
    //      PDE_.SetBCs();
    //
    //      algsys_->ConstructEffectiveMatrix(matrix_factor_);
    //
    //      algsys_->BuildInDirichlet();
    //
    //      // put mass and damping on RHS
    //      TS_alg_->UpdateRHS(actSol);
    //
    //      algsys_->RemoveIDBCInfoFromMatrix();
    //
    //      // substract K^* u^k from RHS
    //      TS_alg_->SubstractStiffnessFromRHS(actSol);
    //
    //      algsys_->SetupPrecond;
    //      algsys_->SetupSolver;
    //      algsys_->Solve;
    //
    //      // new solution is only an increment of the full solution =============
    //      algsys_->GetSolutionVal( solInc );
    //      Double residualL2Norm;
    //      Double etaLineSearch = 1.0;
    //
    //      residualL2Norm = solInc.NormL2();
    //
    //      if ( lineSearch_ == "none" ) {
    //        actSol.Add( 1.0, solInc );
    //      }
    //      else {
    //        residualL2Norm = LineSearchMaterial(solInc, actSol, etaLineSearch, RhsLinL2Norm);
    //      }
    //
    //      residualL2Norm = solInc.NormL2();
    //
    //      PDE_.SaveSolution( actSol );
    //
    //      SBM_Vector actRHS;
    //      algsys_->GetRHSVal( actRHS );
    //
    //      Vector<Double> u_uOld(uOld.GetSize());
    //      u_uOld.Init();
    //      for (UInt ii=0;ii<uOld.GetSize();ii++){
    //        if(uOld[ii]!=0)
    //          u_uOld[ii]=(actSol[ii]-uOld[ii])/uOld[ii];
    //
    //      }
    //      Double incrementL2Norm = u_uOld.NormL2();
    //      std::cout<<"-- residual2Norm = " << residualL2Norm
    //               <<", incrementL2Norm = "<<incrementL2Norm<< std::endl;
    //
    //      Double residualErr;
    //      if ( RhsLinL2Norm > 1.0 )
    //        residualErr    = residualL2Norm /  RhsLinL2Norm;
    //      else
    //        residualErr    = residualL2Norm;
    //
    //      // calculate incremental error
    //      Double solIncrL2Norm = solInc.NormL2();
    //      Double actSolL2Norm = actSol.NormL2();
    //      Double incrementalErr;
    //
    //      if ( actSolL2Norm > 1.0)
    //        incrementalErr = solIncrL2Norm / actSolL2Norm;
    //      else
    //        incrementalErr = solIncrL2Norm;
    //
    //
    //      // --------------------------------------------------------------------
    //      // output of norms and data
    //      // --------------------------------------------------------------------
    //      if ( nonLinLogging_ == true )
    //        WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);
    //
    //      // boolean variable, holds condition if another iteration step
    //      // is necessary
    //      performOneMoreStep =
    //        (incrementL2Norm > incStopCrit_)||(residualErr > residualStopCrit_);
    //
    //    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
    //
    //    // perform corrector step
    //    TS_alg_->Corrector(actSol);
    
  }
  
  
  
  void StdSolveStep::PostStepTrans( ) {
    
    //    WARN("Biot-Savart not yet included";)
    //    Vector<Double> & solHelp =
    //      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    //
    //    // Following method is essential for fractional damping model
    //    TS_alg_->AdvanceTimestep(solHelp);
    //    
    //    // check for Biot Savart
    //    if ( PDE_.IsBiotSavart() ) {
    //      Vector<Double> & sol = 
    //          dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    //      Vector<Double>& magVecBiotSavart = 
    //          PDE_.GetBiotSavart()->CalcFieldAllEqns(false);
    //      sol += magVecBiotSavart;
    //    }
  }
  
  
  
  // ======================================================
  // Solve Step Harmonic  SECTION
  // ======================================================w
  
  void StdSolveStep::PreStepHarmonic() {
    algsys_->InitRHS();
  }
  
  
  void StdSolveStep::SolveStepHarmonic() {
    if ( nonLin_ ) {
      StepHarmonicNonLin();
    }
    else {
      StepHarmonicLin();
    }
  }
  
  
  void StdSolveStep::StepHarmonicLin() {
    
    //JUST A HACK!!!!
    //matrix_factor_Complex_[NO_FCT_ID][STIFFNESS] = Complex(1.0,0);
    //matrix_factor_Complex_[NO_FCT_ID][DAMPING] = CompSetLinRHSlex(0.0,actFreq_*2*M_PI);
    //matrix_factor_Complex_[NO_FCT_ID][MASS] = Complex(-1.0 * actFreq_*actFreq_*4*M_PI*M_PI,0);
    
    //matrix_factor_Complex_[NO_FCT_ID][STIFFNESS] = Complex(1.0,0);
    //matrix_factor_Complex_[NO_FCT_ID][DAMPING] = Complex(1.0,0.0);
    //matrix_factor_Complex_[NO_FCT_ID][MASS] = Complex(1.0,0.0);
    
    
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
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }
    
    algsys_->Solve();
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
  
  void StdSolveStep::GetEigenMode( UInt numMode ) {
    
    
    algsys_->GetEigenMode( numMode );
    
    // Get the solution and store it
    algsys_->GetSolutionVal(solVec_);
  }
  
  
  
  
  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================
  
  // sets excitation coil and returns L2Norm of them
  Double StdSolveStep::SetLinRHS( Double loadFactor, bool nonlin)
  {
    
    //std::cout << "SetLinRHS with bool nonlin = " << nonlin << std::endl;
    Double RhsLinL2Norm;
    
    // to incorporate loads+
    if(nonlin){
      assemble_->AssembleNonLinRHS();
    } else {
      assemble_->AssembleLinRHS();
    }
    
    //Set special RHS Values
    PDE_.SetRhsValues();
    
    // Stores rhs vector into extForces and returns that L2-norm
    algsys_->GetRHSVal( RhsLinVal_ );
    
    RhsLinVal_.ScalarMult(loadFactor);
    
    RhsLinL2Norm = RhsLinVal_.NormL2();
    
    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      // Note: there are PDEs, such as elecconduction, which always have rhs=0. Those should not emit this warning. SE.
      if (pdename_ != "elecConduction")
        WARN("Zero external force vector!! ");
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
    const Double eta[nrEtas] = {0.1, 0.25, 0.5, 1.0}; //, 0.5, 0.25, 0.125, 0.1};
    
    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;
    
    for( UInt i=0; i<nrEtas; i++) {
      //std::cout << "Testing eta = " << eta[i] << std::endl;
      actSol.Add( 1.0, solOld, eta[i], solIncrement);
      
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
    }
    
    //std::cout << "Optimal eta = " << etaOpt << std::endl;
    etaLineSearch = etaOpt;
    
    // Set new solution
    actSol.Add( 1.0, solOld, etaOpt, solIncrement );
    
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
  void StdSolveStep::ReadNonLinData() {
    
    //		std::cout << "Read Non Lin Data" << std::endl;
    // Get ParamNode of pde
    PtrParamNode nonLinNode = solStrat_->GetNonLinNode();
    
    // Check, if any nonlinear node was found
    if( !nonLinNode ) {
      WARN("Taking default parameters for nonlinear data" );
    }
    
    // Read data, if "nonLinear" element was found
    if( nonLinNode ) {
      
      // solution method
      nonLinNode->GetValue( "method", nonLinMethod_, ParamNode::PASS );
      
      // perform logging?
      nonLinNode->GetValue( "logging", nonLinLogging_, ParamNode::PASS );
      
      // type of line search
      if( nonLinNode->Has("lineSearch") ) {
        nonLinNode->Get( "lineSearch")
        ->GetValue( "type", lineSearch_,ParamNode::PASS );
      }
      
      // incremental stopping criterion
      nonLinNode->GetValue( "incStopCrit", incStopCrit_, ParamNode::PASS );
      
      // residual stopping criterion
      nonLinNode->GetValue( "resStopCrit", residualStopCrit_, 
              ParamNode::PASS );
      
      // for hysteresis
      UInt evalParameter;
      nonLinNode->GetValue( "hysteresisEvaluationParameter", evalParameter, ParamNode::INSERT );
      
      /*
       * New: evaluationParameter is a multidigit integer that is used to set
       * multiple parameter at once;
       * it is decomposed into blocks of one or more digits which then act
       * as parameter
       *
       * evalParameter = zyx...dcba
       *  a: evalVersion
       *      0: fixPoint for debugging (hyst has no effect on solution)
       *      1: fixPoint
       *      2: fixPoint with extended rhs
       *      3: deltaMat
       *      4: deltaMat with extended rhs
       *  b: forceReevaluation
       *      0: incremental evaluation of residual
       *      1: reevaluation of residual
       *  c: debugFlag
       *      0: no output
       *      1: debug output to cout
       *
       *  zyx...d > get passed to hystoperator where it is decomposed further
       *
       */
      //      std::cout << "Read in hysteresisEvaluationParameter: " << evalParameter << std::endl;
      evalVersion_ = evalParameter%10;
      
      evalParameter = evalParameter/10;
      forceReevaluation_ = bool(evalParameter%10);
      
      evalParameter = evalParameter/10;
      debugOutput_ = bool(evalParameter%10);
      
      remainingEvalParameter_ = evalParameter/10;
      
      //      std::cout << "evalVersion: " << evalVersion_ << std::endl;
      //      std::cout << "forceReevaluation: " << forceReevaluation_ << std::endl;
      //      std::cout << "debugOutput: " << debugOutput_ << std::endl;
      //      std::cout << "remainingEvalParameter: " << remainingEvalParameter_ << std::endl;
      
      // set parameter for coefFncHyst
      //			std::cout << "set flag in coeffnchyst" << std::endl;
      PDE_.SetFlagInCoefFncHyst("SetInputDependentFlags",remainingEvalParameter_);
      
      //TODO: add defaultValue to xmlSchema!!!!!
      
      //      nonLinNode->GetValue( "evalVersion", evalVersion_, ParamNode::INSERT );
      //   //   std::cout << "evalVersion: " << evalVersion_ << std::endl;
      //      nonLinNode->GetValue( "forceResidualReevaluation", forceReevaluation_, ParamNode::INSERT );
      //   //   std::cout << "forceReevaluation_: " << forceReevaluation_ << std::endl;
      
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
        solutionString = nonLinNode->Get( "stopOnLimit")
                ->Get( "quantity" )->As<std::string>();
        solutionLimit_ = SolutionTypeEnum.Parse(solutionString);
        if (feFunctions_.find(solutionLimit_) == feFunctions_.end() ) 
          EXCEPTION("ERROR: Solution type '" << solutionString << "' is not part of PDE '" << pdename_ << 
                  "' and cannot serve as stopping limit criterion");
        // minimum value
        nonLinNode->Get( "stopOnLimit")
        ->GetValue( "min", minValidValue_, ParamNode::PASS );
        // maximum value
        nonLinNode->Get( "stopOnLimit")
        ->GetValue( "max", maxValidValue_, ParamNode::PASS );
        // region
        std::string solutionRegion;
        solutionRegion = nonLinNode->Get( "stopOnLimit")
                ->Get( "region" )->As<std::string>();
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
  }
  
  void StdSolveStep::WriteNonLinIterToInfoXML(const std::string& pdeName,
          const UInt solStep,
          const UInt iterationCounter,
          const Double residualErr, 
          const Double incrementalErr, double etaLineSearch)
  {
    
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    iter = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT)
            ->Get("iteration",ParamNode::APPEND);
    iter->Get("pdeName")->SetValue(pdeName);
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("residualErr")->SetValue(residualErr);
    iter->Get("incrementalErr")->SetValue(incrementalErr);
    if(etaLineSearch)
      iter->Get("eta_linesearch")->SetValue(etaLineSearch);
  }
  
  void StdSolveStep::WriteNonLinIterToInfoXML(const std::string& pdeName,
          const UInt coupledIterStep,
          const UInt solStep,
          const UInt iterationCounter,
          const Double residualErr, 
          const Double incrementalErr, double etaLineSearch)
  {
    
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    iter = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT);
    iter = iter->GetByVal("couplingStep", "value", coupledIterStep, ParamNode::INSERT)
            ->Get("iteration",ParamNode::APPEND);
    iter->Get("pdeName")->SetValue(pdeName);
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("residualErr")->SetValue(residualErr);
    iter->Get("incrementalErr")->SetValue(incrementalErr);
    // SE: include PDE name and time step
    if(etaLineSearch)
      iter->Get("eta_linesearch")->SetValue(etaLineSearch);
  }
  
} // end of namespace
