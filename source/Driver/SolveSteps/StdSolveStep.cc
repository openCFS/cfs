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
#include "Domain/Results/MHTimeFreqResult.hh"


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
    // Since the entries of solVec_ are pointers to the SingleVector
    // of the FE function, it automatically inserts the values there
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
      EXCEPTION("Time stepping for hysteresis no longer implemented in stdsolvestep");
    }
    
    //currently not supported
    //  else if ( nonLin_ && nonLinMaterial_ ) {
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
      
     // Since the entries of solVec_ are pointers to the SingleVector
     // of the FE function, it automatically inserts the values there
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
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }
    
    algsys_->Solve();
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

    // Set some variables
    UInt N = solStrat_->GetNumHarmN();
    UInt M = solStrat_->GetNumHarmM();
    Double bF = solStrat_->GetBaseFreq();
    UInt numFFT = solStrat_->GetNumFFT();
    bool zeroHarm = solStrat_->IsZeroHarm();
    if(numFFT % 2 != 0){
      EXCEPTION("Please provide a numFFT xml attribute, which is even!");
    }

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
    */
    //std::map<FEMatrixType,Double> empty;
    //algsys_->ConstructEffectiveMatrix(NO_FCT_ID,  empty, true );


    // Incorporate Boundary conditions and
    // recalculate the preconditioner eventually
    algsys_->BuildInDirichlet();

    algsys_->SetupPrecond();
    algsys_->SetupSolver();

    // Solve the linear system
    algsys_->Solve();


    // Get the solution of the initial (linear) multiharmonic system.
    solVecMH_.ResetEntryType(BaseMatrix::EntryType::COMPLEX);
    algsys_->GetFullMultiHarmSolutionVal( solVecMH_, false);

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
    do {
      iterationCounter++;
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
      algsys_->SetupPrecond();
      algsys_->SetupSolver();

      // Solve the deflect system K(u^k) \cdot \Delta u^{k+1} = f - K(u^k) \cdot u^k
      // for the deflect-vector \Delta u^{k+1}
      // TODO DO WE NEED TO CALL IT WITH SETIDBC?
      algsys_->Solve();
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
        if (PDE_.IsIterCoupled()) {
          WriteNonLinIterToInfoXML(pdename_, couplingIter_, 1, iterationCounter, residualErr, incrementalErr, etaLineSearch);
        } else {
          WriteNonLinIterToInfoXML(pdename_, 1, iterationCounter, residualErr, incrementalErr, etaLineSearch);
        }
        // write norm to file
        logFile_ <<  iterationCounter << "\t"
            << residualErr << "\t"
            << incrementalErr << "\t"
            << etaLineSearch << std::endl;
      }

      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
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

    std::cout << "  - Calculating BiLinearForms for multiharmonic analysis" <<std::endl;

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

        // Ok, now it gets confusing because in the performance-optimized
        // version, we draw a border between the harmonics of the system matrix
        // and the solution vector
        Integer tmpH = -solStrat_->GetNumHarmN() + 2 * i;
        Integer h = (tmpH < 0)? tmpH - 1 : tmpH + 1;
        mParser_->SetValue(MathParser::GLOB_HANDLER, "harmonicHandle", h);

        // And set the corresponding frequency
        Double tmpf = (Double)solStrat_->GetBaseFreq();
        Double freq = tmpf * h;
        mParser_->SetValue(MathParser::GLOB_HANDLER, "f", freq);

        if (IS_LOG_ENABLED(stdsolvestep, dbg3)) {
          std::cout<<"harmonic = "<<h<<", frequency = "<<freq<<std::endl;
        }
        if (std::abs(h) > (Integer)M || h == 0) {
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

  void StdSolveStep::GetEigenMode( UInt numMode ) {
    
    
    algsys_->GetEigenMode( numMode );
    
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



    const UInt nrEtas = 4;
    const Double eta[nrEtas] = {0.1, 0.25, 0.5, 1.0}; //, 0.5, 0.25, 0.125, 0.1};

    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    for( UInt i=0; i<nrEtas; i++) {
      LOG_DBG(stdsolvestep) <<" LineSearchMultHarm: Testing eta = " << eta[i];

      if(actSol.GetEntryType() == BaseMatrix::DOUBLE){
        EXCEPTION("StdSolveStep::LineSearchMultHarm Solution vector is real valued in multiharmonic analysis!");
      }else{
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

    }
    
    //std::cout << "Optimal eta = " << etaOpt << std::endl;
    etaLineSearch = etaOpt;

    // Set new solution
    actSol.Add( (Complex)1.0, solOld, etaOpt, solIncrement );
    
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
    
    //    std::cout << "Read Non Lin Data" << std::endl;
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
      
      //TODO: add defaultValue to xmlSchema!!!!!
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
